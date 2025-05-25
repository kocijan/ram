// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <deque>
#include <stdexcept>
#include <iostream>

namespace ram {

MinimizerEngine::MinimizerEngine(
    std::shared_ptr<thread_pool::ThreadPool> thread_pool,
    std::uint32_t k,
    std::uint32_t w,
    std::uint32_t bandwidth,
    std::uint32_t chain,
    std::uint32_t matches,
    std::uint32_t gap)
    : k_(std::min(std::max(k, 1U), 63U)),
      w_(w),
      bandwidth_(bandwidth),
      chain_(chain),
      matches_(matches),
      gap_(gap),
      occurrence_(-1),
      index_(1U << std::min(14U, 2 * k_)),
      thread_pool_(thread_pool ?
          thread_pool :
          std::make_shared<thread_pool::ThreadPool>(1)) {}

std::uint32_t MinimizerEngine::Index::Find(
    std::uint64_t key,
    const Kmer** dst) const {
  auto it = locator.find(key << 1);
  if (it == locator.end()) {
    return 0;
  }
  if (it->first & 1) {
    *dst = &(it->second);
    return 1;
  }
  *dst = &(kmers[it->second.origin >> 32]);
  return static_cast<std::uint32_t>(it->second.origin);
}

void MinimizerEngine::Minimize(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    bool minhash,
    bool hpc) {

  for (auto& it : index_) {
    it.kmers.clear();
    it.locator.clear();
  }

  if (first >= last) {
    return;
  }

  std::vector<std::vector<Kmer>> minimizers(index_.size());
  {
    std::uint64_t mask = index_.size() - 1;

    while (first != last) {
      std::size_t batch_size = 0;
      std::vector<std::future<std::vector<Kmer>>> futures;
      for (; first != last && batch_size < 50000000; ++first) {
        batch_size += (*first)->inflated_len;
        futures.emplace_back(thread_pool_->Submit(
            [&] (decltype(first) it) -> std::vector<Kmer> {
              return Minimize(*it, minhash, hpc);
            },
            first));
      }
      for (auto& it : futures) {
        for (const auto& jt : it.get()) {
          auto& m = minimizers[jt.value() & mask];
          if (m.capacity() == m.size()) {
            m.reserve(m.capacity() * 1.5);
          }
          m.emplace_back(jt);
        }
      }
    }
  }

  {
    std::vector<std::future<std::pair<std::size_t, std::size_t>>> futures;
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      futures.emplace_back(thread_pool_->Submit(
          [&] (std::uint32_t i) -> std::pair<std::size_t, std::size_t> {
            if (minimizers[i].empty()) {
              return std::make_pair(0, 0);
            }

            RadixSort(
                minimizers[i].begin(),
                minimizers[i].end(),
                62U,
                Kmer::SortByValue);

                        // stop dummy
            minimizers[i].emplace_back(~minimizers[i].back().description, -1);

            std::size_t num_origins = 0;
            std::size_t num_keys = 0;

            for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {  // NOLINT
              if (minimizers[i][j - 1].value() != minimizers[i][j].value()) {
                if (c > 1) {
                  num_origins += c;
                }
                ++num_keys;
                c = 0;
              }
            }

            return std::make_pair(num_origins, num_keys);
          },
          i));
    }
    for (std::uint32_t i = 0; i < minimizers.size(); ++i) {
      auto num_entries = futures[i].get();
      if (minimizers[i].empty()) {
        continue;
      }

      index_[i].kmers.reserve(num_entries.first);
      index_[i].locator.reserve(num_entries.second);

      for (std::uint64_t j = 1, c = 1; j < minimizers[i].size(); ++j, ++c) {
        if (minimizers[i][j - 1].value() != minimizers[i][j].value()) {
          if (c == 1) {
            index_[i].locator.emplace(
                minimizers[i][j - 1].value() << 1 | 1,
                minimizers[i][j - 1]);
          } else {
            index_[i].locator.emplace(
                minimizers[i][j - 1].value() << 1,
                Kmer(0, index_[i].kmers.size() << 32 | c));
            for (std::uint64_t k = j - c; k < j; ++k) {
              index_[i].kmers.emplace_back(minimizers[i][k]);
            }
          }
          c = 0;
        }
      }

      std::vector<Kmer>().swap(minimizers[i]);
    }
  }
}

void MinimizerEngine::Filter(double frequency) {
  if (!(0 <= frequency && frequency <= 1)) {
    throw std::invalid_argument(
        "[ram::MinimizerEngine::Filter] error: invalid frequency");
  }

  if (frequency == 0) {
    occurrence_ = -1;
    return;
  }

  std::vector<std::uint32_t> occurrences;
  for (const auto& it : index_) {
    for (const auto& jt : it.locator) {
      if (jt.first & 1) {
        occurrences.emplace_back(1);
      } else {
        occurrences.emplace_back(static_cast<std::uint32_t>(jt.second.origin));
      }
    }
  }

  if (occurrences.empty()) {
    occurrence_ = -1;
    return;
  }

  std::nth_element(
      occurrences.begin(),
      occurrences.begin() + (1 - frequency) * occurrences.size(),
      occurrences.end());
  occurrence_ = occurrences[(1 - frequency) * occurrences.size()] + 1;
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    bool avoid_equal,
    bool avoid_symmetric,
    bool minhash,
    bool hpc,
    std::vector<std::uint32_t>* filtered) const {
  auto sketch = Minimize(sequence, minhash, hpc);
  if (sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  std::vector<Match> matches;
  auto add_match = [&] (const Kmer& kmer, const Kmer* origin) -> void {

    if (avoid_equal && sequence->id == origin->id()) {
      return;
    }
    if (avoid_symmetric && sequence->id > origin->id()) {
      return;
    }

    std::uint64_t rhs_id = origin->id();
    std::uint64_t strand_ = kmer.strand() == origin->strand();
    std::uint64_t lhs_pos = kmer.position();
    std::uint16_t lhs_span = kmer.span();
    std::uint64_t rhs_pos = origin->position();
    std::uint16_t rhs_span = origin->span();
    std::uint64_t diagonal = !strand_ ?
        rhs_pos + lhs_pos :
        rhs_pos - lhs_pos + (3ULL << 30);

    // Print for match: the id of the sequence, the start and end positions of the match in the sequence, the id of the other sequence, and the start and end positions of the match in the other sequence, the length of the match in the sequence, the length of the match in the other sequence, the diagonal of the match, and the strand of the match.
    // std::cerr << "match " << sequence->id << " " << lhs_pos << " " << lhs_pos + lhs_span << " " << rhs_id << " " << (strand_ ? rhs_pos : rhs_pos + rhs_span) << " " << (strand_ ? rhs_pos + rhs_span : rhs_pos) << " " << lhs_span << " " << rhs_span << " " << diagonal << " " << strand_ << std::endl;

    matches.emplace_back(
        (((rhs_id << 1) | strand_) << 32) | diagonal,
        (lhs_pos << 32) | rhs_pos,
        (lhs_span << 8) | rhs_span);
  };

  struct Hit {
    const Kmer* kmer;
    std::uint32_t n;
    const Kmer* origins;

    Hit(const Kmer* kmer, std::uint32_t n, const Kmer* origins)
        : kmer(kmer),
          n(n),
          origins(origins) {}

    bool operator<(const Hit& other) const {
      return n < other.n;
    }
  };
  std::vector<Hit> filtered_hits;

  std::uint64_t mask = index_.size() - 1;
  std::uint32_t prev = 0;

  sketch.emplace_back(-1, sequence->inflated_len << 1);  // stop dummy

  for (const auto& kmer : sketch) {
    std::uint32_t i = kmer.value() & mask;
    const Kmer* origins = nullptr;
    auto n = index_[i].Find(kmer.value(), &origins);
    if (n > occurrence_) {
      filtered_hits.emplace_back(&kmer, n, origins);
      if (filtered) {
        filtered->emplace_back(kmer.position());
      }
      continue;
    }

    std::size_t rescuees = std::min(
        static_cast<std::size_t>(kmer.position() - prev) / bandwidth_,
        filtered_hits.size());
    if (rescuees) {
      std::partial_sort(
          filtered_hits.begin(),
          filtered_hits.begin() + rescuees,
          filtered_hits.end());
      for (auto it = filtered_hits.begin(); rescuees; rescuees--, ++it) {
        for (; it->n; it->n--, ++it->origins) {
          add_match(*it->kmer, it->origins);
        }
      }
    }
    filtered_hits.clear();
    prev = kmer.position();

    for (; n; n--, ++origins) {
      add_match(kmer, origins);
    }
  }

  // Group matches by (rhs_id, strand) pairs and call ChainDP separately
  std::vector<biosoup::Overlap> all_overlaps;

  if (!matches.empty())
  {
    RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);

    std::uint64_t current_group = matches[0].group >> 32;
    std::size_t group_start = 0;

    for (std::size_t i = 1; i <= matches.size(); ++i)
    {
      if (i == matches.size() || (matches[i].group >> 32) != current_group)
      {
        // Extract matches for current group
        std::vector<Match> group_matches;
        group_matches.reserve(i - group_start);
        for (std::size_t j = group_start; j < i; ++j)
        {
          group_matches.emplace_back(std::move(matches[j]));
        }

        // Call ChainDP for this group
        auto group_overlaps = ChainDP(sequence->id, std::move(group_matches));

        // Move results to combined vector
        all_overlaps.insert(all_overlaps.end(),
                            std::make_move_iterator(group_overlaps.begin()),
                            std::make_move_iterator(group_overlaps.end()));

        if (i < matches.size())
        {
          current_group = matches[i].group >> 32;
          group_start = i;
        }
      }
    }
  }

  return all_overlaps;
}

std::vector<biosoup::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid>& lhs,
    const std::unique_ptr<biosoup::NucleicAcid>& rhs,
    bool minhash,
    bool hpc) const {

  auto lhs_sketch = Minimize(lhs, minhash, hpc);
  if (lhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  auto rhs_sketch = Minimize(rhs, minhash, hpc);
  if (rhs_sketch.empty()) {
    return std::vector<biosoup::Overlap>{};
  }

  RadixSort(lhs_sketch.begin(), lhs_sketch.end(), 62U, Kmer::SortByValue);
  RadixSort(rhs_sketch.begin(), rhs_sketch.end(), 62U, Kmer::SortByValue);

  std::uint64_t rhs_id = rhs->id;

  std::vector<Match> matches;
  for (std::uint32_t i = 0, j = 0; i < lhs_sketch.size(); ++i) {
    while (j < rhs_sketch.size()) {
      if (lhs_sketch[i].value() < rhs_sketch[j].value()) {
        break;
      } else if (lhs_sketch[i].value() == rhs_sketch[j].value()) {
        for (std::uint32_t k = j; k < rhs_sketch.size(); ++k) {
          if (lhs_sketch[i].value() != rhs_sketch[k].value()) {
            break;
          }

          std::uint64_t strand =
              (lhs_sketch[i].strand() & 1) == (rhs_sketch[k].strand() & 1);
          std::uint64_t lhs_pos = lhs_sketch[i].position();
          std::uint16_t lhs_span = lhs_sketch[i].span();
          std::uint64_t rhs_pos = rhs_sketch[k].position();
          std::uint16_t rhs_span = rhs_sketch[k].span();
          std::uint64_t diagonal = !strand ?
              rhs_pos + lhs_pos :
              rhs_pos - lhs_pos + (3ULL << 30);

          matches.emplace_back(
              (((rhs_id << 1) | strand) << 32) | diagonal,
              (lhs_pos << 32) | rhs_pos,
              (lhs_span << 8) | rhs_span);
        }
        break;
      } else {
        ++j;
      }
    }
  }

  // Group matches by (rhs_id, strand) pairs and call ChainDP separately
  std::vector<biosoup::Overlap> all_overlaps;

  if (!matches.empty())
  {
    RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);

    std::uint64_t current_group = matches[0].group >> 32;
    std::size_t group_start = 0;

    for (std::size_t i = 1; i <= matches.size(); ++i)
    {
      if (i == matches.size() || (matches[i].group >> 32) != current_group)
      {
        // Extract matches for current group
        std::vector<Match> group_matches;
        group_matches.reserve(i - group_start);
        for (std::size_t j = group_start; j < i; ++j)
        {
          group_matches.emplace_back(std::move(matches[j]));
        }

        // Call ChainDP for this group
        auto group_overlaps = ChainDP(lhs->id, std::move(group_matches));

        // Move results to combined vector
        all_overlaps.insert(all_overlaps.end(),
                            std::make_move_iterator(group_overlaps.begin()),
                            std::make_move_iterator(group_overlaps.end()));

        if (i < matches.size())
        {
          current_group = matches[i].group >> 32;
          group_start = i;
        }
      }
    }
  }

  return all_overlaps;
}

std::vector<biosoup::Overlap> MinimizerEngine::Chain(
    std::uint64_t lhs_id,
    std::vector<Match>&& matches) const {
  RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);
  matches.emplace_back(-1, -1, -1);  // stop dummy

  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].group - matches[j].group > bandwidth_) {
      // if (i - j >= 4) {
      if (i - j >= std::min(chain_, 4U))
      {

        if (!intervals.empty() && intervals.back().second > j)
        { // extend
          intervals.back().second = i;
        }
        else
        { // new
          intervals.emplace_back(j, i);
        }
      }
      ++j;
      while (j < i && matches[i].group - matches[j].group > bandwidth_) {
        ++j;
      }
    }
  }

  std::vector<biosoup::Overlap> dst;
  for (const auto& it : intervals) {
    std::uint64_t j = it.first;
    std::uint64_t i = it.second;

    if (i - j < chain_) {
      continue;
    }

    RadixSort(
        matches.begin() + j,
        matches.begin() + i,
        64,
        Match::SortByPositions);

    std::uint64_t strand = matches[j].strand();

    std::vector<std::uint64_t> indices;
    if (strand) {  // same strand
      indices = LongestSubsequence(  // increasing
          matches.begin() + j,
          matches.begin() + i,
          std::less<std::uint64_t>());
    } else {  // different strand
      indices = LongestSubsequence(  // decreasing
          matches.begin() + j,
          matches.begin() + i,
          std::greater<std::uint64_t>());
    }

    if (indices.size() < chain_) {
      continue;
    }

    indices.emplace_back(matches.size() - 1 - j);  // stop dummy from above
    for (std::uint64_t k = 1, l = 0; k < indices.size(); ++k) {
      if (matches[j + indices[k]].lhs_position() -
          matches[j + indices[k - 1]].lhs_position() > gap_) {
        if (k - l < chain_) {
          l = k;
          continue;
        }

        std::uint32_t lhs_matches = 0;
        std::uint32_t lhs_begin = 0;
        std::uint32_t lhs_end = 0;
        std::uint32_t rhs_matches = 0;
        std::uint32_t rhs_begin = 0;
        std::uint32_t rhs_end = 0;

        for (std::uint64_t m = l; m < k; ++m) {
          const auto& match = matches[j + indices[m]];
          std::uint32_t lhs_pos = match.lhs_position();
          if (lhs_pos > lhs_end) {
            lhs_matches += lhs_end - lhs_begin;
            lhs_begin = lhs_pos;
          }
          lhs_end = lhs_pos + match.lhs_span();

          std::uint32_t rhs_pos = match.rhs_position();
          rhs_pos = strand ?
              rhs_pos :
              (1U << 31) - (rhs_pos + match.rhs_span() - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + match.rhs_span();
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < matches_) {
          l = k;
          continue;
        }

        /// Print for chain: the id of the sequence, the start and end positions of the chain in the sequence, the id of the other sequence, and the start and end positions of the chain in the other sequence, the length of the chain in the sequence, the length of the chain in the other sequence, and the strand of the chain.
        // std::cerr << "chain " << lhs_id << " " << matches[j + indices[l]].lhs_position() << " " << matches[j + indices[k - 1]].lhs_position() + matches[j + indices[k - 1]].lhs_span() << " " << matches[j + indices[l]].rhs_id() << " " << (strand ? matches[j + indices[l]].rhs_position() : matches[j + indices[k - 1]].rhs_position()) << " " << (strand ? matches[j + indices[k - 1]].rhs_position() + matches[j + indices[k - 1]].rhs_span() : matches[j + indices[l]].rhs_position() + matches[j + indices[l]].rhs_span()) << " " << std::min(lhs_matches, rhs_matches) << " " << strand << std::endl;

        const auto& first = matches[j + indices[l]];
        const auto& last  = matches[j + indices[k - 1]];

        dst.emplace_back(
            lhs_id,
            first.lhs_position(),
             last.lhs_position() + last.lhs_span(),
            matches[j].rhs_id(),
            strand ?
                first.rhs_position() :
                 last.rhs_position(),
            strand ?
                 last.rhs_position() +  last.rhs_span() :
                first.rhs_position() + first.rhs_span(),
            std::min(lhs_matches, rhs_matches),
            strand);

        l = k;
      }
    }
  }
  return dst;
}

static inline float log2f(float x) // NB: this doesn't work when x<2
{
  union
  {
    float f;
    uint32_t i;
  } z = {x};
  float log_2 = ((z.i >> 23) & 255) - 128;
  z.i &= ~(255 << 23);
  z.i += 127 << 23;
  log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
  return log_2;
}

// Helper function to compute score between two anchors
int32_t MinimizerEngine::ComputeDPScore(
    const std::pair<std::uint32_t, std::uint32_t> &ai,
    const std::pair<std::uint32_t, std::uint32_t> &aj,
    std::uint32_t max_dist_x,
    std::uint32_t max_dist_y,
    std::uint32_t bandwidth,
    float chain_gap_scale,
    float chain_skip_scale) const
{

  int32_t dq = ai.second - aj.second; // Distance in query
  int32_t dr = ai.first - aj.first;   // Distance in reference

  if (dq <= 0 || dq > static_cast<int32_t>(max_dist_x))
  {
    return std::numeric_limits<int32_t>::min();
  }

  if (dr <= 0 || dq > static_cast<int32_t>(max_dist_y))
  {
    return std::numeric_limits<int32_t>::min();
  }

  int32_t dd = std::abs(dr - dq); // Deviation from diagonal
  if (dd > static_cast<int32_t>(bandwidth))
  {
    return std::numeric_limits<int32_t>::min();
  }

  int32_t dg = std::min(dr, dq); // Min of distances
  int32_t q_span = k_;           // Use the kmer length

  int32_t score = std::min(q_span, dg); // Base score

  if (dd || dg > q_span)
  {
    float lin_pen = chain_gap_scale * dd + chain_skip_scale * dg;
    float log_pen = dd >= 1 ? log2f(dd + 1) : 0.0f;
    score -= static_cast<int32_t>(lin_pen + 0.5f * log_pen);
  }

  return score;
}

// New DP-based chaining function with same signature as Chain
std::vector<biosoup::Overlap> MinimizerEngine::ChainDP(
    std::uint64_t lhs_id,
    std::vector<Match> &&matches) const
{
  // Set default parameters
  std::uint32_t max_dist_x = 500;          // Max gap in reference
  std::uint32_t max_dist_y = 500;          // Max gap in query
  std::uint32_t dp_bandwidth = bandwidth_; // Use existing bandwidth parameter
  std::uint32_t max_skip = 25;             // Max anchors to skip
  std::uint32_t max_iter = 5000;           // Max iterations
  int32_t min_cnt = chain_;                // Min anchors in chain, use existing parameter
  int32_t min_sc = matches_;               // Min score required, use existing parameter
  float chain_gap_scale = 0.008f * k_;     // Gap cost scale
  float chain_skip_scale = 0.000f * k_;    // Skip cost scale
  int32_t max_drop = dp_bandwidth;         // Max score drop

  if (matches.empty())
  {
    return std::vector<biosoup::Overlap>{};
  }

  if (max_dist_x < dp_bandwidth)
  {
    max_dist_x = dp_bandwidth;
  }
  if (max_dist_y < dp_bandwidth)
  {
    max_dist_y = dp_bandwidth;
  }

  // Sort matches by group (to identify strands)
  RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);
  matches.emplace_back(-1, -1, -1); // stop dummy

  std::vector<biosoup::Overlap> overlaps;

  // Skip interval identification and process all matches together
  // Check if we have enough matches to form a chain
  if (matches.size() - 1 < static_cast<std::size_t>(min_cnt))
  {
    return overlaps;
  }

  // Sort matches by position
  RadixSort(
      matches.begin(),
      matches.end() - 1, // Exclude the dummy
      64,
      Match::SortByPositions);

  std::uint64_t strand = matches[0].strand();

  // Convert matches to anchors for DP
  std::vector<std::pair<std::uint32_t, std::uint32_t>> anchors;
  for (std::size_t k = 0; k < matches.size() - 1; ++k)
  {
    // Only include matches with the same strand orientation
    if (matches[k].strand() == strand)
    {
      anchors.emplace_back(
          matches[k].rhs_position(),
          matches[k].lhs_position());
    }
  }

  int64_t n_a = anchors.size();
  if (n_a < min_cnt)
  {
    return overlaps;
  }

  // Allocate memory for DP arrays
  std::vector<int32_t> f(n_a);     // Score array
  std::vector<int64_t> p(n_a, -1); // Predecessor array
  std::vector<int32_t> t(n_a, 0);  // Temporary array for backtracking

  // Fill the score and backtrack arrays
  for (int64_t a_i = 0, st = 0; a_i < n_a; ++a_i)
  {
    int64_t max_j = -1;
    int32_t max_f = 1; // Default score for a single anchor
    uint32_t n_skip = 0;

    // Find appropriate starting point
    while (st < a_i && (anchors[a_i].first - anchors[st].first > max_dist_x))
    {
      ++st;
    }

    // Limit iterations
    if (a_i - st > static_cast<int64_t>(max_iter))
    {
      st = a_i - max_iter;
    }

    // DP calculation - find best predecessor
    for (int64_t a_j = a_i - 1; a_j >= st; --a_j)
    {
      int32_t sc = ComputeDPScore(
          anchors[a_i], anchors[a_j],
          max_dist_x, max_dist_y,
          dp_bandwidth, chain_gap_scale,
          chain_skip_scale);

      if (sc == std::numeric_limits<int32_t>::min())
      {
        continue;
      }

      sc += f[a_j];
      if (sc > max_f)
      {
        max_f = sc;
        max_j = a_j;
        if (n_skip > 0)
        {
          --n_skip;
        }
      }
      else if (t[a_j] == a_i)
      {
        if (++n_skip > max_skip)
        {
          break;
        }
      }

      if (p[a_j] >= 0)
      {
        t[p[a_j]] = a_i;
      }
    }

    // Set score and predecessor
    f[a_i] = max_f;
    p[a_i] = max_j;
  }

  // Backtrack to find chains - process by score (highest to lowest)
  std::vector<bool> used(n_a, false);
  std::vector<std::vector<int64_t>> chains;

  // Sort anchors by score in descending order
  std::vector<std::pair<int32_t, int64_t>> scored_anchors;
  for (int64_t i = 0; i < n_a; ++i)
  {
    if (f[i] >= min_sc)
    {
      scored_anchors.emplace_back(f[i], i);
    }
  }

  std::sort(scored_anchors.begin(), scored_anchors.end(),
            std::greater<std::pair<int32_t, int64_t>>());

  // Process anchors in decreasing order of score
  for (const auto &scored_anchor : scored_anchors)
  {
    int64_t a_i = scored_anchor.second;

    if (used[a_i])
      continue;

    std::vector<int64_t> chain;
    int32_t max_drop_so_far = 0;
    int32_t max_f = f[a_i];
    int64_t curr = a_i;

    // Backtrack to build the chain
    while (curr >= 0)
    {
      if (used[curr])
        break;

      used[curr] = true;
      chain.push_back(curr);

      if (p[curr] >= 0)
      {
        int32_t gap_sc = f[curr] - f[p[curr]];
        max_drop_so_far = std::max(max_drop_so_far, gap_sc);
        if (max_drop_so_far > max_drop)
          break;
      }

      curr = p[curr];
    }

    // Add chain if it meets criteria
    if (chain.size() >= static_cast<size_t>(min_cnt) && max_f >= min_sc)
    {
      std::reverse(chain.begin(), chain.end());
      chains.push_back(chain);
    }
  }

  // Convert chains to overlaps
  for (const auto &chain : chains)
  {
    if (chain.empty())
      continue;

    std::uint32_t lhs_matches = 0;
    std::uint32_t lhs_begin = 0;
    std::uint32_t lhs_end = 0;
    std::uint32_t rhs_matches = 0;
    std::uint32_t rhs_begin = 0;
    std::uint32_t rhs_end = 0;

    // Calculate match lengths
    for (const auto &idx : chain)
    {
      const auto &match_idx = matches[idx]; // Get actual match from the original array

      std::uint32_t lhs_pos = match_idx.lhs_position();
      if (lhs_pos > lhs_end)
      {
        lhs_matches += lhs_end - lhs_begin;
        lhs_begin = lhs_pos;
      }
      lhs_end = lhs_pos + match_idx.lhs_span();

      std::uint32_t rhs_pos = match_idx.rhs_position();
      rhs_pos = strand ? rhs_pos : (1U << 31) - (rhs_pos + match_idx.rhs_span() - 1);
      if (rhs_pos > rhs_end)
      {
        rhs_matches += rhs_end - rhs_begin;
        rhs_begin = rhs_pos;
      }
      rhs_end = rhs_pos + match_idx.rhs_span();
    }
    lhs_matches += lhs_end - lhs_begin;
    rhs_matches += rhs_end - rhs_begin;

    if (std::min(lhs_matches, rhs_matches) < matches_)
    {
      continue;
    }

    // Create the overlap
    const auto &first_match = matches[chain.front()];
    const auto &last_match = matches[chain.back()];

    // Print chain details
    // std::cerr << "dp_chain " << lhs_id << " "
    //           << first_match.lhs_position() << " "
    //           << last_match.lhs_position() + last_match.lhs_span() << " "
    //           << first_match.rhs_id() << " "
    //           << (strand ? first_match.rhs_position() : last_match.rhs_position()) << " "
    //           << (strand ? last_match.rhs_position() + last_match.rhs_span() : first_match.rhs_position() + first_match.rhs_span()) << " "
    //           << std::min(lhs_matches, rhs_matches) << " "
    //           << strand << std::endl;

    overlaps.emplace_back(
        lhs_id,
        first_match.lhs_position(),
        last_match.lhs_position() + last_match.lhs_span(),
        first_match.rhs_id(),
        strand ? first_match.rhs_position() : last_match.rhs_position(),
        strand ? last_match.rhs_position() + last_match.rhs_span() : first_match.rhs_position() + first_match.rhs_span(),
        std::min(lhs_matches, rhs_matches),
        strand);

    // break?
    break;
  }

  return overlaps;
}

std::vector<MinimizerEngine::Kmer> MinimizerEngine::Minimize(
    const std::unique_ptr<biosoup::NucleicAcid>& sequence,
    bool minhash,
    bool hpc) const {
  if (sequence->inflated_len < k_) {
    return std::vector<Kmer>{};
  }

  std::uint64_t mask = (1ULL << k_) - 1;

  auto hash = [&] (std::uint64_t key) -> std::uint64_t {
    key = ~key + (key << 21);
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key = key ^ key >> 28;
    key = key + (key << 31);
    return key;
  };

  std::deque<Kmer> window;
  auto window_add = [&] (std::uint64_t description, std::uint64_t origin) -> void {
    while (!window.empty() && window.back().value() > (description >> 8)) {
      window.pop_back();
    }
    window.emplace_back(description, origin);
  };
  auto window_update = [&] (std::uint32_t position) -> void {
    while (!window.empty() && window.front().position() < position) {
      window.pop_front();
    }
  };

  std::uint64_t shift = k_ - 1;
  std::uint64_t minimizer_lo = 0;
  std::uint64_t minimizer_hi = 0;
  std::uint64_t reverse_minimizer_lo = 0;
  std::uint64_t reverse_minimizer_hi = 0;
  std::uint64_t id = static_cast<std::uint64_t>(sequence->id) << 32;
  std::uint64_t is_stored = 1ULL;

  std::vector<Kmer> dst;

  for (std::uint32_t i = 0, win_span = 0, kmer_span = 0, base_cnt = 0; i < sequence->inflated_len; ++i, ++win_span, ++kmer_span) {
    std::uint64_t c = sequence->Code(i);
    
    // skip homopolymer
    if (hpc && i && sequence->Code(i - 1) == c) {
      continue;
    }
    // found new char
    base_cnt++;

    // remove last from kmer
    if (base_cnt > k_) {
      kmer_span--;
      if (hpc) {
        auto last_c = sequence->Code(i - kmer_span - 1);
        while (sequence->Code(i - kmer_span) == last_c) kmer_span--;
      }
    }

    minimizer_lo = ((minimizer_lo << 1) | (c & 1)) & mask;
    minimizer_hi = ((minimizer_hi << 1) | (c & 2)) & mask;
    reverse_minimizer_lo = (reverse_minimizer_lo >> 1) | (((c ^ 3) & 1) << shift);
    reverse_minimizer_hi = (reverse_minimizer_hi >> 1) | (((c ^ 3) & 2) << shift);
    if (base_cnt >= k_ && kmer_span < 256ULL) {
      std::uint64_t origin =
      (std::uint64_t(i + 1U - kmer_span) << 33) |
      ((base_cnt - (k_ - 1U)) << 1);
      if (minimizer_hi < reverse_minimizer_hi) {
        window_add(((hash(minimizer_lo) + hash(minimizer_hi)) << 8) | kmer_span, origin);
      } else if (minimizer_hi > reverse_minimizer_hi) {
        origin |= 1ULL << 32;
        window_add(((hash(reverse_minimizer_lo) + hash(reverse_minimizer_hi)) << 8) | kmer_span, origin);
      }
    }
    if (base_cnt >= (k_) + (w_ - 1U)) {
      for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->value() != window.front().value()) {
          break;
        }
        if (it->origin & is_stored) {
          continue;
        }
        dst.emplace_back(it->description, id | (it->origin >> 32));
        it->origin |= is_stored;
      }
      win_span--;
      if (hpc) {
        auto last_c = sequence->Code(i - win_span - 1);
        while (sequence->Code(i - win_span) == last_c) win_span--;
      }
      window_update(base_cnt - 1U - (k_ - 1U) - (w_ - 1U) + 1U);
    }
  }

  if (minhash) {
    RadixSort(dst.begin(), dst.end(), 62U, Kmer::SortByValue);
    dst.resize(sequence->inflated_len / k_);
    RadixSort(dst.begin(), dst.end(), 64U, Kmer::SortByOrigin);
  }

  return dst;
}

template<typename RandomAccessIterator, typename Compare>
void MinimizerEngine::RadixSort(
    RandomAccessIterator first,
    RandomAccessIterator last,
    std::uint8_t max_bits,
    Compare comp) {  //  unary comparison function
  if (first >= last) {
    return;
  }

  std::vector<typename std::iterator_traits<RandomAccessIterator>::value_type> tmp(last - first);  // NOLINT
  auto begin = tmp.begin();
  auto end = tmp.end();

  std::uint64_t buckets[0x100]{};  // 256 b
  std::uint8_t shift = 0;
  for (; shift < max_bits; shift += 8) {
    std::uint64_t counts[0x100]{};
    for (auto it = first; it != last; ++it) {
      ++counts[comp(*it) >> shift & 0xFF];
    }
    for (std::uint64_t i = 0, j = 0; i < 0x100; j += counts[i++]) {
      buckets[i] = j;
    }
    for (auto it = first; it != last; ++it) {
      *(begin + buckets[comp(*it) >> shift & 0xFF]++) = *it;
    }
    std::swap(begin, first);
    std::swap(end, last);
  }

  if (shift / 8 & 1) {  // copy the sorted array for odd cases
    for (; first != last; ++first, ++begin) {
      *begin = *first;
    }
  }
}

template<typename Compare>
std::vector<std::uint64_t> MinimizerEngine::LongestSubsequence(
    std::vector<Match>::const_iterator first,
    std::vector<Match>::const_iterator last,
    Compare comp) {  // binary comparison function
  if (first >= last) {
    return std::vector<std::uint64_t>{};
  }

  std::vector<std::uint64_t> minimal(last - first + 1, 0);
  std::vector<std::uint64_t> predecessor(last - first, 0);

  std::uint64_t longest = 0;
  for (auto it = first; it != last; ++it) {
    std::uint64_t lo = 1, hi = longest;
    while (lo <= hi) {
      std::uint64_t mid = lo + (hi - lo) / 2;
      if ((first + minimal[mid])->lhs_position() < it->lhs_position() &&
          comp((first + minimal[mid])->rhs_position(), it->rhs_position())) {
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    predecessor[it - first] = minimal[lo - 1];
    minimal[lo] = it - first;
    longest = std::max(longest, lo);
  }

  std::vector<std::uint64_t> dst;
  for (std::uint64_t i = 0, j = minimal[longest]; i < longest; ++i) {
    dst.emplace_back(j);
    j = predecessor[j];
  }
  std::reverse(dst.begin(), dst.end());

  return dst;
}

}  // namespace ram
