// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <deque>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "biosoup/timer.hpp"
#include "biosoup/progress_bar.hpp"

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
        kmer_frequency_cutoff_(0),
        index_(1U << std::min(14U, 2 * k_)),
        thread_pool_(thread_pool ? thread_pool : std::make_shared<thread_pool::ThreadPool>(1)) {}

  std::uint32_t MinimizerEngine::Index::Find(
      std::uint64_t key,
      const Kmer **dst) const
  {
    auto it = locator.find(key << 1);
    if (it == locator.end())
    {
      return 0;
    }
    if (it->first & 1)
    {
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

  biosoup::Timer timer{};

  // Step 1: Count k-mers in all sequences
  std::cerr << "[ram::MinimizerEngine::Minimize] counting k-mers..." << std::endl;
  timer.Start();

  kmer_frequencies_ = CountKmers(first, last, hpc);

  std::cerr << "[ram::MinimizerEngine::Minimize] counted " << kmer_frequencies_.size()
            << " unique k-mers " << std::fixed << timer.Stop() << "s" << std::endl;

  // Step 2: Estimate cutoff for frequent k-mers (use same logic as existing Filter)
  timer.Start();
  std::cerr << "[ram::MinimizerEngine::Minimize] calculating frequency cutoff..." << std::endl;

  if (!kmer_frequencies_.empty())
  {
    std::vector<std::uint32_t> counts;
    counts.reserve(kmer_frequencies_.size());
    for (const auto &kv : kmer_frequencies_)
    {
      counts.push_back(kv.second);
    }

    // Use 0.001 as default frequency threshold for k-mer filtering
    double kmer_frequency_threshold = 0.01;
    std::nth_element(
        counts.begin(),
        counts.begin() + (1 - kmer_frequency_threshold) * counts.size(),
        counts.end());
    kmer_frequency_cutoff_ = counts[(1 - kmer_frequency_threshold) * counts.size()] + 1;
  }

  std::cerr << "[ram::MinimizerEngine::Minimize] frequency cutoff: " << kmer_frequency_cutoff_
            << " " << std::fixed << timer.Stop() << "s" << std::endl;

  // Step 3: Generate minimizers with frequency filtering
  timer.Start();
  std::cerr << "[ram::MinimizerEngine::Minimize] generating minimizers..." << std::endl;

  std::vector<std::vector<Kmer>> minimizers(index_.size());
  {
    std::uint64_t mask = index_.size() - 1;
    auto first_copy = first;

    while (first_copy != last)
    {
      std::size_t batch_size = 0;
      std::vector<std::future<std::vector<Kmer>>> futures;
      for (; first_copy != last && batch_size < 50000000; ++first_copy)
      {
        batch_size += (*first_copy)->inflated_len;
        futures.emplace_back(thread_pool_->Submit(
            [&](decltype(first_copy) it) -> std::vector<Kmer>
            {
              return Minimize(*it, minhash, hpc, &kmer_frequencies_, kmer_frequency_cutoff_);
            },
            first_copy));
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

  std::cerr << "[ram::MinimizerEngine::Minimize] generated minimizers "
            << std::fixed << timer.Stop() << "s" << std::endl;

  // Step 4: Build index
  timer.Start();
  std::cerr << "[ram::MinimizerEngine::Minimize] building index..." << std::endl;

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

  std::cerr << "[ram::MinimizerEngine::Minimize] built index "
            << std::fixed << timer.Stop() << "s" << std::endl;
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

std::vector<ram::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid> &sequence,
    bool avoid_equal,
    bool avoid_symmetric,
    bool minhash,
    bool hpc,
    std::vector<std::uint32_t> *filtered) const
{
  auto sketch = Minimize(sequence, minhash, hpc);
  if (sketch.empty()) {
    return std::vector<ram::Overlap>{};
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

    matches.emplace_back(
        (((rhs_id << 1) | strand_) << 32) | diagonal,
        (lhs_pos << 32) | rhs_pos,
        (lhs_span << 8) | rhs_span,
        kmer.value()); // Store k-mer hash
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
        static_cast<std::size_t>(kmer.position() - prev) / 500,
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
  std::vector<ram::Overlap> all_overlaps;

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

std::vector<ram::Overlap> MinimizerEngine::Map(
    const std::unique_ptr<biosoup::NucleicAcid> &lhs,
    const std::unique_ptr<biosoup::NucleicAcid> &rhs,
    bool minhash,
    bool hpc) const
{

  auto lhs_sketch = Minimize(lhs, minhash, hpc);
  if (lhs_sketch.empty()) {
    return std::vector<ram::Overlap>{};
  }

  auto rhs_sketch = Minimize(rhs, minhash, hpc);
  if (rhs_sketch.empty()) {
    return std::vector<ram::Overlap>{};
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
              (lhs_span << 8) | rhs_span,
              lhs_sketch[i].value()); // Store k-mer hash from lhs
        }
        break;
      } else {
        ++j;
      }
    }
  }

  // Group matches by (rhs_id, strand) pairs and call ChainDP separately
  std::vector<ram::Overlap> all_overlaps;

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

std::vector<ram::Overlap> MinimizerEngine::Chain(
    std::uint64_t lhs_id,
    std::vector<Match> &&matches) const
{
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

  std::vector<ram::Overlap> dst;
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
            strand ? first.rhs_position() : last.rhs_position(),
            strand ? last.rhs_position() + last.rhs_span() : first.rhs_position() + first.rhs_span(),
            std::min(lhs_matches, rhs_matches),
            strand,
            0); // dp_score = 0 for LIS-based chaining

        l = k;
      }
    }
  }
  return dst;
}

// New DP-based chaining function with same signature as Chain
std::vector<ram::Overlap> MinimizerEngine::ChainDP(
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
  int32_t max_drop = dp_bandwidth;         // Max score drop
  float bandwidth_threshold = 0.05f;       // Threshold for bandwidth
  float bandwidth_penalty = 1.0 / bandwidth_threshold; // Penalty for exceeding bandwidth

  if (matches.empty())
  {
    return std::vector<ram::Overlap>{};
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

  std::vector<ram::Overlap> overlaps;

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
  matches.emplace_back(-1, -1, -1); // stop dummy

  std::uint64_t strand = matches[0].strand();

  // Convert matches to anchors for DP
  std::vector<std::pair<std::uint32_t, std::uint32_t>> anchors;
  std::vector<std::int32_t> anchor_frequencies;
  for (std::size_t k = 0; k < matches.size() - 1; ++k)
  {
    anchors.emplace_back(
        matches[k].rhs_position(),
        matches[k].lhs_position());

    // Extract actual k-mer frequency from the stored k-mer hash
    std::int32_t frequency = 1; // Default frequency

    if (!kmer_frequencies_.empty())
    {
      auto it = kmer_frequencies_.find(matches[k].kmer_hash);
      if (it != kmer_frequencies_.end())
      {
        frequency = static_cast<std::int32_t>(it->second);
      }
    }

    anchor_frequencies.emplace_back(frequency);
  }

  int64_t n_a = anchors.size();
  if (n_a < min_cnt)
  {
    return overlaps;
  }

  // Allocate memory for DP arrays
  std::vector<int32_t> f(n_a);     // Score array
  std::vector<int64_t> p(n_a, -1); // Predecessor array
  std::vector<int64_t> indels(n_a, 0);      // Cumulative indel counts
  std::vector<int64_t> self_length(n_a, 0); // Cumulative self lengths
  std::vector<int32_t> t(n_a, -1);          // Temporary array for backtracking

  // Fill the score and backtrack arrays
  for (int64_t a_i = 0, st = 0; a_i < n_a; ++a_i)
  {
    int64_t max_j = -1;
    uint32_t n_chn_skip = 0;
    uint32_t n_max_skip = 0;
    int64_t max_indels = 0;
    int64_t max_self_length = 0;
    int32_t max_score = (int32_t)k_ >= anchor_frequencies[a_i] ? (int32_t)k_ / anchor_frequencies[a_i] : 1; // Normalize by frequency

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
      int64_t distance_pos = anchors[a_i].first - anchors[a_j].first;
      int64_t distance_self_pos = strand ? anchors[a_i].second - anchors[a_j].second : anchors[a_j].second - anchors[a_i].second;
      // Skip invalid transitions
      if (distance_pos > max_dist_x || distance_self_pos > max_dist_y)
      {
        continue;
      }
      if (distance_pos <= 0 || distance_self_pos <= 0)
      {
        continue;
      }
      // Calculate gap and cumulative values
      int64_t distance_gap = std::abs(distance_self_pos - distance_pos);
      int64_t total_indels = indels[a_j] + distance_gap;
      int64_t total_self_length = self_length[a_j] + distance_self_pos;
      // Check bandwidth constraints
      if (total_indels > bandwidth_threshold * total_self_length)
      {
        continue;
      }
      // Calculate score
      int64_t distance_min = std::min(distance_self_pos, distance_pos);
      int32_t score = std::min(
          (int64_t)k_,
          distance_min);
      score = score >= anchor_frequencies[a_j] ? score / anchor_frequencies[a_j] : 1; // Normalize by frequency

      // Apply gap rate penalty
      float gap_rate = total_indels / static_cast<float>(total_self_length);
      score -= (int32_t)(gap_rate * score * bandwidth_penalty);

      score += f[a_j]; // Add predecessor score

      if (score > max_score)
      {
        max_score = score;
        max_j = a_j;
        max_indels = total_indels;
        max_self_length = total_self_length;
        n_max_skip = 0;
        if (n_chn_skip > 0)
        {
          --n_chn_skip;
        }
      }
      else
      {
        if (++n_max_skip > max_skip)
        {
          break;
        }
        if (t[a_j] == a_i)
        {
          if (++n_chn_skip > max_skip)
          {
            break;
          }
        }
      }

      if (p[a_j] >= 0)
      {
        t[p[a_j]] = a_i;
      }
    }

    // Set score and predecessor
    f[a_i] = max_score;
    p[a_i] = max_j;
    indels[a_i] = max_indels;
    self_length[a_i] = max_self_length;
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

      if (p[curr] >= 0 && !used[p[curr]])
      {
        for (int64_t k = p[curr] + 1; k < curr; ++k)
        {
          used[k] = true;
        }
      }

      curr = p[curr];
    }

    // Add chain if it meets criteria
    if (chain.size() >= static_cast<size_t>(min_cnt) && max_f >= min_sc)
    {
      if (strand)
        std::reverse(chain.begin(), chain.end());
      chains.push_back(chain);
    }
  }

  // Convert chains to overlaps
  std::uint32_t prev_matches = 0;
  for (std::size_t chain_idx = 0; chain_idx < chains.size(); ++chain_idx)
  {
    const auto &chain = chains[chain_idx];
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

    std::uint32_t current_matches = std::min(lhs_matches, rhs_matches);

    // Break if this is not the first iteration and matches are less than 90% of previous
    if (prev_matches != 0 && current_matches < static_cast<std::uint32_t>(0.9 * prev_matches))
    {
      break;
    }

    // Create the overlap
    const auto &first_match = matches[chain.front()];
    const auto &last_match = matches[chain.back()];

    // Get dp_score from the last element of the chain (which is the first when we built it backwards)
    int32_t chain_dp_score = strand ? f[chain.back()] : f[chain.front()];

    overlaps.emplace_back(
        lhs_id,
        first_match.lhs_position(),
        last_match.lhs_position() + last_match.lhs_span(),
        first_match.rhs_id(),
        strand ? first_match.rhs_position() : last_match.rhs_position(),
        strand ? last_match.rhs_position() + last_match.rhs_span() : first_match.rhs_position() + first_match.rhs_span(),
        std::min(lhs_matches, rhs_matches),
        strand,
        chain_dp_score);

    prev_matches = current_matches;
  }

  return overlaps;
}

MinimizerEngine::KmerFrequencyTable MinimizerEngine::CountKmers(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
    std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
    bool hpc) const
{
  std::uint64_t total_sequences = std::distance(first, last);

  // Process only 10sqrt(n) reads for faster k-mer frequency estimation
  std::uint64_t sample_size = static_cast<std::uint64_t>(std::sqrt(total_sequences) * 10.);
  sample_size = std::max(sample_size, static_cast<std::uint64_t>(1)); // Ensure at least 1 sequence
  sample_size = std::min(sample_size, total_sequences);               // Don't exceed total

  std::cerr << "[ram::MinimizerEngine::CountKmers] processing " << sample_size
            << " out of " << total_sequences << " sequences (sqrt sampling)" << std::endl;

  // Pre-allocate with estimated size to reduce rehashing
  KmerFrequencyTable frequency_table;
  frequency_table.reserve(sample_size * 1000); // Reduced estimate based on sample size

  std::uint64_t processed_sequences = 0;
  std::uint64_t next_milestone = sample_size / 10;
  if (next_milestone == 0)
    next_milestone = 1;

  biosoup::Timer timer{};
  timer.Start();

  // Calculate step size for even distribution
  std::uint64_t step = total_sequences / sample_size;
  if (step == 0)
    step = 1;

  // Use fixed batch size based on sample count rather than total count
  const std::uint64_t batch_size = std::max(static_cast<std::uint64_t>(1), sample_size / (thread_pool_->num_threads() * 4));

  auto current_it = first;
  std::uint64_t sequences_processed = 0;

  while (sequences_processed < sample_size && current_it != last)
  {
    std::vector<std::future<KmerFrequencyTable>> futures;

    // Create batches with sampled sequences
    for (std::uint32_t i = 0; i < thread_pool_->num_threads() && sequences_processed < sample_size; ++i)
    {
      std::vector<std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator> batch_sequences;
      std::uint64_t sequences_in_batch = std::min(batch_size, sample_size - sequences_processed);

      // Collect sequences for this batch using step sampling
      for (std::uint64_t j = 0; j < sequences_in_batch && current_it != last; ++j)
      {
        batch_sequences.push_back(current_it);

        // Advance by step size for even distribution
        for (std::uint64_t k = 0; k < step && current_it != last; ++k)
        {
          ++current_it;
        }
      }

      if (!batch_sequences.empty())
      {
        futures.emplace_back(thread_pool_->Submit(
            [this, hpc](std::vector<std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator> batch_seqs) -> KmerFrequencyTable
            {
              KmerFrequencyTable local_table;
              // Pre-allocate based on expected k-mers per sequence
              local_table.reserve(batch_seqs.size() * 1000);

              for (auto seq_it : batch_seqs)
              {
                CountKmersInSequence(*seq_it, local_table, hpc);
              }
              return local_table;
            },
            batch_sequences));

        sequences_processed += batch_sequences.size();
      }
    }

    // Collect results and merge sequentially to avoid contention
    std::vector<KmerFrequencyTable> local_results;
    local_results.reserve(futures.size());

    for (auto &future : futures)
    {
      local_results.emplace_back(future.get());
    }

    // Sequential merge to avoid hash table contention
    for (const auto &local_table : local_results)
    {
      for (const auto &kv : local_table)
      {
        frequency_table[kv.first] += kv.second;
      }
    }

    processed_sequences = sequences_processed;

    // Progress reporting
    if (processed_sequences >= next_milestone || processed_sequences >= sample_size)
    {
      double progress = (double)processed_sequences / sample_size * 100.0;
      std::cerr << "[ram::MinimizerEngine::CountKmers] processed " << processed_sequences
                << "/" << sample_size << " sequences (" << std::fixed << std::setprecision(1)
                << progress << "%) " << std::fixed << std::setprecision(3) << timer.Lap() << "s" << std::endl;

      next_milestone += sample_size / 10;
      if (next_milestone > sample_size)
      {
        next_milestone = sample_size;
      }
    }
  }

  std::cerr << "[ram::MinimizerEngine::CountKmers] completed k-mer counting on " << processed_sequences
            << " sequences " << std::fixed << timer.Stop() << "s" << std::endl;

  // Print histogram of k-mer frequencies
  if (!frequency_table.empty())
  {
    std::cerr << "\n[ram::MinimizerEngine::CountKmers] K-mer Frequency Distribution:" << std::endl;
    std::cerr << "============================================================" << std::endl;

    // Count how many k-mers have each frequency (1 to 100)
    std::vector<std::uint64_t> histogram(101, 0); // Index 0 unused, 1-100 for frequencies
    std::uint64_t max_frequency = 0;
    std::uint64_t total_kmers = 0;

    for (const auto &kv : frequency_table)
    {
      std::uint32_t freq = kv.second;
      max_frequency = std::max(max_frequency, static_cast<std::uint64_t>(freq));
      total_kmers++;

      if (freq <= 100)
      {
        histogram[freq]++;
      }
    }

    std::cerr << "Total unique k-mers: " << total_kmers << std::endl;
    std::cerr << "Maximum frequency: " << max_frequency << std::endl;
    std::cerr << "\nFrequency Distribution (top 100):" << std::endl;
    std::cerr << "Freq | Count    | Percentage | Histogram" << std::endl;
    std::cerr << "-----|----------|------------|--------------------------------------------------" << std::endl;

    // Find max count for scaling the histogram bars
    std::uint64_t max_count = *std::max_element(histogram.begin() + 1, histogram.end());
    const int bar_width = 40;

    for (int freq = 1; freq <= 100; ++freq)
    {
      if (histogram[freq] > 0)
      {
        double percentage = (double)histogram[freq] / total_kmers * 100.0;
        int bar_length = max_count > 0 ? (int)((double)histogram[freq] / max_count * bar_width) : 0;

        std::cerr << std::setw(4) << freq << " | "
                  << std::setw(8) << histogram[freq] << " | "
                  << std::setw(9) << std::fixed << std::setprecision(3) << percentage << "% | ";

        // Draw histogram bar
        for (int i = 0; i < bar_length; ++i)
        {
          std::cerr << "█";
        }
        std::cerr << std::endl;
      }
    }

    // Count k-mers with frequency > 100
    std::uint64_t high_freq_count = 0;
    for (const auto &kv : frequency_table)
    {
      if (kv.second > 100)
      {
        high_freq_count++;
      }
    }

    if (high_freq_count > 0)
    {
      double percentage = (double)high_freq_count / total_kmers * 100.0;
      std::cerr << ">100 | " << std::setw(8) << high_freq_count << " | "
                << std::setw(9) << std::fixed << std::setprecision(3) << percentage << "% | "
                << "(frequencies > 100)" << std::endl;
    }

    std::cerr << "============================================================" << std::endl;

    // Print some summary statistics
    std::uint64_t singletons = histogram[1];
    std::uint64_t doubletons = histogram[2];
    double singleton_rate = (double)singletons / total_kmers * 100.0;
    double doubleton_rate = (double)doubletons / total_kmers * 100.0;

    std::cerr << "\nSummary Statistics:" << std::endl;
    std::cerr << "- Singletons (freq=1): " << singletons << " ("
              << std::fixed << std::setprecision(1) << singleton_rate << "%)" << std::endl;
    std::cerr << "- Doubletons (freq=2): " << doubletons << " ("
              << std::fixed << std::setprecision(1) << doubleton_rate << "%)" << std::endl;
    std::cerr << "- High frequency (>100): " << high_freq_count << " ("
              << std::fixed << std::setprecision(1) << (double)high_freq_count / total_kmers * 100.0 << "%)" << std::endl;
    std::cerr << std::endl;
  }

  return frequency_table;
}

void MinimizerEngine::CountKmersInSequence(
    const std::unique_ptr<biosoup::NucleicAcid> &sequence,
    KmerFrequencyTable &frequency_table,
    bool hpc) const
{
  if (sequence->inflated_len < k_)
  {
    return;
  }

  std::uint64_t mask = (1ULL << k_) - 1;

  // Use same hash function as in Minimize
  auto hash = [&](std::uint64_t key) -> std::uint64_t
  {
    key = ~key + (key << 21);
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key = key ^ key >> 28;
    key = key + (key << 31);
    return key;
  };

  std::uint64_t minimizer_lo = 0;
  std::uint64_t minimizer_hi = 0;
  std::uint64_t reverse_minimizer_lo = 0;
  std::uint64_t reverse_minimizer_hi = 0;
  std::uint64_t shift = k_ - 1;

  // Pre-increment frequency_table load factor to reduce rehashing
  frequency_table.reserve(frequency_table.size() + sequence->inflated_len / k_);

  for (std::uint32_t i = 0, kmer_span = 0, base_cnt = 0; i < sequence->inflated_len; ++i, ++kmer_span)
  {
    std::uint64_t c = sequence->Code(i);

    // skip homopolymer
    if (hpc && i && sequence->Code(i - 1) == c)
    {
      continue;
    }
    // found new char
    base_cnt++;

    // remove last from kmer
    if (base_cnt > k_)
    {
      kmer_span--;
      if (hpc)
      {
        auto last_c = sequence->Code(i - kmer_span - 1);
        while (sequence->Code(i - kmer_span) == last_c)
          kmer_span--;
      }
    }

    minimizer_lo = ((minimizer_lo << 1) | (c & 1)) & mask;
    minimizer_hi = ((minimizer_hi << 1) | (c & 2)) & mask;
    reverse_minimizer_lo = (reverse_minimizer_lo >> 1) | (((c ^ 3) & 1) << shift);
    reverse_minimizer_hi = (reverse_minimizer_hi >> 1) | (((c ^ 3) & 2) << shift);

    if (base_cnt >= k_ && kmer_span < 256ULL)
    {
      // Count k-mer using same hash combination as in Minimize
      std::uint64_t kmer_hash;
      if (minimizer_hi < reverse_minimizer_hi)
      {
        kmer_hash = hash(minimizer_lo) + hash(minimizer_hi);
      }
      else if (minimizer_hi > reverse_minimizer_hi)
      {
        kmer_hash = hash(reverse_minimizer_lo) + hash(reverse_minimizer_hi);
      }
      else
      {
        continue; // Skip palindromic k-mers
      }

      ++frequency_table[kmer_hash]; // Use pre-increment for slight performance gain
    }
  }
}

std::vector<MinimizerEngine::Kmer> MinimizerEngine::Minimize(
    const std::unique_ptr<biosoup::NucleicAcid> &sequence,
    bool minhash,
    bool hpc,
    const KmerFrequencyTable *frequency_table,
    std::uint32_t frequency_cutoff) const
{

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
  auto window_add = [&](std::uint64_t description, std::uint64_t origin) -> void
  {
    // Check frequency filter if provided
    if (frequency_table && frequency_cutoff > 0)
    {
      std::uint64_t kmer_hash = description >> 8;
      auto it = frequency_table->find(kmer_hash);
      if (it != frequency_table->end() && it->second >= frequency_cutoff)
      {
        return; // Skip this k-mer as it's too frequent
      }
    }

    // Get frequency for new k-mer
    std::uint64_t new_kmer_hash = description >> 8;
    std::uint32_t new_frequency = 1;
    if (frequency_table)
    {
      auto it = frequency_table->find(new_kmer_hash);
      if (it != frequency_table->end())
      {
        new_frequency = it->second;
      }
    }

    // Compare with window back based on frequency (keep less frequent k-mers)
    while (!window.empty())
    {
      std::uint64_t window_kmer_hash = window.back().value();
      std::uint32_t window_frequency = 1;
      if (frequency_table)
      {
        auto it = frequency_table->find(window_kmer_hash);
        if (it != frequency_table->end())
        {
          window_frequency = it->second;
        }
      }

      // If new k-mer is less frequent (lower frequency count), remove window back
      // If frequencies are equal, fall back to hash comparison
      if (new_frequency < window_frequency ||
          (new_frequency == window_frequency && new_kmer_hash < window_kmer_hash))
      {
        window.pop_back();
      }
      else
      {
        break;
      }
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
