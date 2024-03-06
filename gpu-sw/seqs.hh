#ifndef __ECCL_SEQS_HH__
#define __ECCL_SEQS_HH__

#include <array>
#include <cstdint>
#include <istream>
#include <cassert>
#include <vector>
#include <memory>
#include <string_view>

#include "core.hh"

namespace eccl {

/*! API for a batch of sequences
 *
 * optimized for medium-sized batches of reads;
 * may be used for a single sequence or a whole genome.
 *
 * limits:
 *   - maximum length of one sequence: 4Gi;
 */


template<seq_type Type>
struct seq_view {
	/*! if reverse/complement sequence is not relevent,
	 *   use bits directly.
	 * otherwise,
	 *   bits&1 means the reverse sequence,
	 *   bits&2 means the complement sequence,
	 *   and use bits without the 2 LSB.
	 */
	const unsigned int* bits;
	/*! sub-seq range: [begin, end).  begin<=end. */
	uint32_t begin, end;
	uint32_t len() const noexcept { return end-begin; }
	std::size_t buf_size() const noexcept {
#ifndef ECCL_COMPACT_CODE
		return Type==seq_type::xna?(len()+7)/8:(len()+3)/4;
#else
		return Type==seq_type::xna?(len()+9)/10:(len()+5)/6;
#endif
	}
};

namespace detail {
struct seqs_base {
	struct Coords {
		uint64_t seq_s;
		uint32_t seq_l;
		uint32_t name_e;
	};
	uint32_t max_seq_l{0};
	std::vector<Coords> coords;
	std::vector<unsigned int> seq_buf;
	std::vector<char> name_buf;
};
template<eccl::seq_type Type, typename Cont>
class BitStream;
}

class fasta_chunk;
template<seq_type Type>
class seqs_builder;

template<seq_type Type>
class seqs {
public:
	explicit seqs(const fasta_chunk& chunk);
	seqs(seqs_builder<Type>&& builder) noexcept;
	~seqs() =default;
	seqs(const seqs&) =delete;
	seqs& operator=(const seqs&) =delete;
	seqs(seqs&& r) noexcept: _b{std::move(r._b)} { }
	seqs& operator=(seqs&& r) noexcept {
		_b=std::move(r._b);
		return *this;
	}
#if 0
	std::pair<uint64_t, uint32_t> start_len(std::size_t i) const noexcept {
		auto& c=_coords[i];
		return {c.seq_s, c.seq_l};
	}
#endif
	std::size_t size() const noexcept { return _b.coords.size(); }
	seq_view<Type> seq(std::size_t idx) const noexcept {
		auto& c=_b.coords[idx];
		return {&_b.seq_buf[c.seq_s], 0, c.seq_l};
	}
	seq_view<Type> operator[](std::size_t idx) const noexcept {
		return seq(idx);
	}
	std::string_view name(std::size_t i) const noexcept {
		auto a=i>0?_b.coords[i-1].name_e:0;
		auto b=_b.coords[i].name_e-a;
		return {&_b.name_buf[a], b};
	}
	uint64_t begin(std::size_t idx) const noexcept {
		auto& c=_b.coords[idx];
		return c.seq_s*code<Type>::n_per_word;
	}

	const unsigned int* buf() const noexcept { return &_b.seq_buf[0]; }
	std::size_t buf_size() const noexcept { return _b.seq_buf.size(); }
	uint32_t max_seq_len() const noexcept { return _b.max_seq_l; }

private:
	detail::seqs_base _b;
	friend std::ostream& operator<<(std::ostream& str, const seqs& seqs) {
		return str<<&seqs._b.seq_buf[0]<<'@'<<seqs._b.coords.size();
	}
};

using xna_seqs=seqs<seq_type::xna>;
using prot_seqs=seqs<seq_type::prot>;

/*! use 3 bits for each DNA/RNA base, but 4-bit aligned in a sequence.
 * use 5 bits for each amino acid.
 */

class seq_encoder_tables;

template<seq_type Type>
class seq_encoder {
public:
	explicit seq_encoder() noexcept =default;
	~seq_encoder() =default;

	code<Type> operator()(char c) const noexcept;
	char rev(code<Type> i) const noexcept;
};

template<seq_type Type>
class seqs_builder {
public:
	seqs_builder(std::size_t num_bases, std::size_t names_size, std::size_t num_seqs);
	~seqs_builder() =default;
	seqs_builder(const seqs_builder&) =delete;
	seqs_builder& operator=(const seqs_builder&) =delete;

	template<bool RC=false>
	void push_back(std::string_view seq);

private:
	detail::seqs_base _b;
	std::shared_ptr<detail::BitStream<Type, std::vector<unsigned int>>> _str;
	seq_encoder<Type> _encoder;
	friend class seqs<Type>;
};

template<seq_type Type>
seqs<Type>::seqs(seqs_builder<Type>&& builder) noexcept:
_b{std::move(builder._b)} { }


class seq_encoder_tables {
public:
	explicit seq_encoder_tables() noexcept;
private:
	std::array<unsigned char, 256> _table_xna, _table_prot;
	std::array<char, 8> _rev_xna;
	std::array<char, 32> _rev_prot;
	template<seq_type Type>
	const unsigned char* table() const noexcept {
		return Type==seq_type::xna?&_table_xna[0]:&_table_prot[0];
	}
	template<seq_type Type>
	const char* rev() const noexcept {
		return Type==seq_type::xna?&_rev_xna[0]:&_rev_prot[0];
	}
	static seq_encoder_tables _t;
	template<seq_type> friend class seq_encoder;
};

template<seq_type Type>
inline code<Type> seq_encoder<Type>::operator()(char c) const noexcept {
	return {seq_encoder_tables::_t.table<Type>()[static_cast<unsigned char>(c)]};
}
template<seq_type Type>
inline char seq_encoder<Type>::rev(code<Type> i) const noexcept {
	return seq_encoder_tables::_t.rev<Type>()[i.value];
}

class fasta_chunk {
public:
	fasta_chunk(fasta_chunk&&) =default;
	~fasta_chunk() =default;
	fasta_chunk(const fasta_chunk&) =delete;
	fasta_chunk& operator=(const fasta_chunk&) =delete;

private:
	explicit fasta_chunk() =default;
	fasta_chunk& operator=(fasta_chunk&&) =default;

	struct Coords {
		std::size_t name_s, name_e;
		std::size_t seq_s, seq_e;
	};
	std::size_t names_total{0};
	std::size_t seqs_total{0};
	std::vector<Coords> coords;
	std::string buf;
	friend class fasta_reader;
	template<seq_type>
	friend class seqs;
};

class fasta_reader {
public:
	explicit fasta_reader() =default;
	~fasta_reader() =default;
	fasta_reader(const fasta_reader&) =delete;
	fasta_reader& operator=(const fasta_reader&) =delete;

	std::shared_ptr<fasta_chunk> load(eccl::source<char>& str, std::size_t buf_size, unsigned int num_seq);

private:
	using Chunk=fasta_chunk;
	Chunk cur{};
	enum class St {
		expect_hdr,
		in_hdr,
		expect_seq,
		in_seq,
		end,
	} st;
	bool hiteof{false};
};

// XXX
class seqs_ref {
public:
private:
};

template<typename Int, typename Alloc, typename ScoreOpts>
void fill_score_table_xna(std::vector<Int, Alloc>& _table, const ScoreOpts& scores) {
	_table.reserve(8*8+2);
	_table.resize(8*8, -100);
	auto set=[&_table](auto a, auto b, Int v) {
		_table[pair<seq_type::xna>{a, b}.value]=v;
		_table[pair<seq_type::xna>{b, a}.value]=v;
	};
	const Int m=scores.match ? *scores.match : 1;
	const Int mm=scores.mismatch ? *scores.mismatch : -4;
	const Int go=scores.gap_open ? *scores.gap_open : -6;
	const Int ge=scores.gap_extend ? *scores.gap_extend : -1;
	_table.push_back(go);
	_table.push_back(ge);

	eccl::seq_encoder<seq_type::xna> b2b{};
	const auto A=b2b('A');
	const auto T=b2b('T');
	const auto C=b2b('C');
	const auto G=b2b('G');
	const auto U=b2b('U');
	set(A, A, m);
	set(T, T, m);
	set(T, U, m);
	set(U, U, m);
	set(C, C, m);
	set(G, G, m);

	set(A, T, mm);
	set(A, U, mm);
	set(A, C, mm);
	set(A, G, mm);

	set(T, C, mm);
	set(T, G, mm);
	set(U, C, mm);
	set(U, G, mm);

	set(C, G, mm);
}
void mat_blosum62(std::array<int, 32*32>& mat);

template<seq_type Type, typename Int, typename Alloc, typename ScoreOpts, typename=std::enable_if_t<Type==seq_type::xna>>
inline void fill_score_table_impl(std::vector<Int, Alloc>& _table, const ScoreOpts& scores) {
	fill_score_table_xna(_table, scores);
}
template<seq_type Type, typename Int, typename Alloc, typename=std::enable_if_t<Type==seq_type::prot>>
inline void fill_score_table_impl(std::vector<Int, Alloc>& _table) {
	_table.reserve(32*32+2);
	std::array<int, 32*32> mat{-100, };
	mat_blosum62(mat);
	for(auto v: mat)
		_table.push_back(v);
	const Int gap_open=-11;
	const Int gap_ext=-1;
	_table.push_back(gap_open);
	_table.push_back(gap_ext);
}
template<seq_type Type, typename Int, typename Alloc, typename ScoreOpts, typename=std::enable_if_t<Type==seq_type::xna>>
inline void fill_score_table(std::vector<Int, Alloc>& _table, const ScoreOpts& scores) {
	return fill_score_table_impl<Type, Int, Alloc>(_table, scores);
}
template<seq_type Type, typename Int, typename Alloc,  typename=std::enable_if_t<Type==seq_type::prot>>
inline void fill_score_table(std::vector<Int, Alloc>& _table) {
	return fill_score_table_impl<Type, Int, Alloc>(_table);
}

}

#endif
