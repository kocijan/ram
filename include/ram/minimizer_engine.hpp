// Copyright (c) 2020 Robet Vaser

#ifndef RAM_MINIMIZER_ENGINE_HPP_
#define RAM_MINIMIZER_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <vector>
#include <unordered_map>
#include <utility>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram {

class MinimizerEngine {
 public:
  MinimizerEngine(
      std::shared_ptr<thread_pool::ThreadPool> thread_pool = nullptr,
      std::uint32_t k = 15,  // element of [1, 63]
      std::uint32_t w = 5,
      std::uint32_t bandwidth = 500,
      std::uint32_t chain = 4,
      std::uint32_t matches = 100,
      std::uint32_t gap = 10000);

  MinimizerEngine(const MinimizerEngine&) = delete;
  MinimizerEngine& operator=(const MinimizerEngine&) = delete;

  MinimizerEngine(MinimizerEngine&&) = default;
  MinimizerEngine& operator=(MinimizerEngine&&) = default;

  ~MinimizerEngine() = default;

  // transform set of sequences to minimizer index
  // minhash = pick only the smallest sequence->data.size() / k minimizers
  void Minimize(
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator first,
      std::vector<std::unique_ptr<biosoup::NucleicAcid>>::const_iterator last,
      bool minhash = false,
      bool hpc = false);

  // set occurrence frequency threshold
  void Filter(double frequency);

  // find overlaps in preconstructed minimizer index
  std::vector<biosoup::Overlap> Map(
      const std::unique_ptr<biosoup::NucleicAcid>& sequence,
      bool avoid_equal,  // ignore overlaps in which lhs_id == rhs_id
      bool avoid_symmetric,  // ignore overlaps in which lhs_id > rhs_id
      bool minhash = false,  // only lhs
      bool hpc = false,  // use homopolymer compression
      std::vector<std::uint32_t>* filtered = nullptr) const;

  // find overlaps between a pair of sequences
  std::vector<biosoup::Overlap> Map(
      const std::unique_ptr<biosoup::NucleicAcid>& lhs,
      const std::unique_ptr<biosoup::NucleicAcid>& rhs,
      bool minhash = false,  // only lhs
      bool hpc = false) const;

 private:
  struct Kmer {
   public:
    Kmer() = default;
    Kmer(std::uint64_t description, std::uint64_t origin)
        : description(description),
          origin(origin) {
    }

    std::uint64_t value() const {
      return description >> 8;
    }

    std::uint8_t span() const {
      return static_cast<std::uint8_t>(description);
    }

    std::uint32_t id() const {
      return static_cast<std::uint32_t>(origin >> 32);
    }

    std::uint32_t position() const {
      return static_cast<std::uint32_t>(origin) >> 1;
    }

    bool strand() const {
      return origin & 1;
    }

    static std::uint64_t SortByValue(const Kmer& kmer) {
      return kmer.value();
    }

    static std::uint64_t SortByOrigin(const Kmer& kmer) {
      return kmer.origin;
    }

    std::uint64_t description;
    std::uint64_t origin;
  };

  struct Match {
   public:
    Match() = default;
    Match(std::uint64_t group, std::uint64_t positions, std::uint16_t spans)
        : group(group),
          positions(positions),
          spans(spans) {
    }

    std::uint32_t rhs_id() const {
      return static_cast<std::uint32_t>(group >> 33);
    }

    bool strand() const {
      return (group >> 32) & 1;
    }

    std::uint32_t diagonal() const {
      return static_cast<std::uint32_t>(group);
    }

    std::uint32_t lhs_position() const {
      return static_cast<std::uint32_t>(positions >> 32);
    }

    std::uint32_t rhs_position() const {
      return static_cast<std::uint32_t>(positions);
    }

    std::uint8_t lhs_span() const {
      return static_cast<std::uint8_t>(spans >> 8);
    }

    std::uint8_t rhs_span() const {
      return static_cast<std::uint8_t>(spans);
    }

    static std::uint64_t SortByGroup(const Match& match) {
      return match.group;
    }
    static std::uint64_t SortByPositions(const Match& match) {
      return match.positions;
    }

    std::uint64_t group;
    std::uint64_t positions;
    std::uint16_t spans;
  };

  class Index {
   public:
    Index() = default;

    std::uint32_t Find(std::uint64_t key, const Kmer** dst) const;

    struct Hash {
      public:
      std::size_t operator()(std::uint64_t key) const {
        return std::hash<std::uint64_t>()(key >> 1);
      }
    };
    struct KeyEqual {
      public:
      bool operator()(std::uint64_t lhs, std::uint64_t rhs) const {
        return (lhs >> 1) == (rhs >> 1);
      }
    };

    std::vector<Kmer> kmers;
    std::unordered_map<std::uint64_t, Kmer, Hash, KeyEqual> locator;
  };

  // hpc = use homopolymer compression
  std::vector<Kmer> Minimize(
      const std::unique_ptr<biosoup::NucleicAcid>& sequence,
      bool minhash = false,
      bool hpc = false) const;

  std::vector<biosoup::Overlap> Chain(
      std::uint64_t lhs_id,
      std::vector<Match>&& matches) const;

  template<typename RandomAccessIterator, typename Compare>
  static void RadixSort(
      RandomAccessIterator first,
      RandomAccessIterator last,
      std::uint8_t max_bits,
      Compare comp);  //  unary comparison function

  template<typename Compare>
  static std::vector<std::uint64_t> LongestSubsequence(
      std::vector<Match>::const_iterator first,
      std::vector<Match>::const_iterator last,
      Compare comp);  // binary comparison function

  std::uint32_t k_;
  std::uint32_t w_;
  std::uint32_t bandwidth_;
  std::uint32_t chain_;
  std::uint32_t matches_;
  std::uint64_t gap_;
  std::uint32_t occurrence_;
  std::vector<Index> index_;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool_;
};

}  // namespace ram

#endif  // RAM_MINIMIZER_ENGINE_HPP_
