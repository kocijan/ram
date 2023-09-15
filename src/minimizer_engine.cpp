// Copyright (c) 2020 Robert Vaser

#include "ram/minimizer_engine.hpp"

#include <deque>
#include <stdexcept>

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

  return Chain(sequence->id, std::move(matches));
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

  return Chain(lhs->id, std::move(matches));
}

std::vector<biosoup::Overlap> MinimizerEngine::Chain(
    std::uint64_t lhs_id,
    std::vector<Match>&& matches) const {
  RadixSort(matches.begin(), matches.end(), 64, Match::SortByGroup);
  matches.emplace_back(-1, -1, -1);  // stop dummy

  std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
  for (std::uint64_t i = 1, j = 0; i < matches.size(); ++i) {  // NOLINT
    if (matches[i].group - matches[j].group > bandwidth_) {
      if (i - j >= 4) {
        if (!intervals.empty() && intervals.back().second > j) {  // extend
          intervals.back().second = i;
        } else {  // new
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
          lhs_end = lhs_pos + k_ + match.lhs_span();

          std::uint32_t rhs_pos = match.rhs_position();
          rhs_pos = strand ?
              rhs_pos :
              (1U << 31) - (rhs_pos + k_ + match.rhs_span() - 1);
          if (rhs_pos > rhs_end) {
            rhs_matches += rhs_end - rhs_begin;
            rhs_begin = rhs_pos;
          }
          rhs_end = rhs_pos + k_+ match.rhs_span();
        }
        lhs_matches += lhs_end - lhs_begin;
        rhs_matches += rhs_end - rhs_begin;
        if (std::min(lhs_matches, rhs_matches) < matches_) {
          l = k;
          continue;
        }

        const auto& first = matches[j + indices[l]];
        const auto& last  = matches[j + indices[k - 1]];

        dst.emplace_back(
            lhs_id,
            first.lhs_position(),
             last.lhs_position() + k_ + last.lhs_span(),
            matches[j].rhs_id(),
            strand ?
                first.rhs_position() :
                 last.rhs_position(),
            strand ?
                 last.rhs_position() + k_ +  last.rhs_span() :
                first.rhs_position() + k_ + first.rhs_span(),
            std::min(lhs_matches, rhs_matches),
            strand);

        l = k;
      }
    }
  }
  return dst;
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
      (std::uint64_t(i - (k_ - 1U) - kmer_span) << 33) |
      ((base_cnt - 1U - (k_ - 1U)) << 1);
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
