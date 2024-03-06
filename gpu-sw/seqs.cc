#include "seqs.hh"

//#include <iostream>
#include <algorithm>
#include <type_traits>


// XXX alternative: 32-bit, 10-bp, 2-bit unused per word
// avoid 4-bit padded to 3-bit conversion.
// but, is 10 harder to divide by than 8???
template<eccl::seq_type Type, typename Cont>
class eccl::detail::BitStream {
public:
	explicit constexpr BitStream(Cont& cont) noexcept:
		cont{cont} { }
	std::size_t size() const noexcept { return tmp_cnt; }
	void add(unsigned int v) {
#ifndef ECCL_COMPACT_CODE
		tmp_bits=(tmp_bits<<padded_width)+v;
		if((++tmp_cnt)%(padded_count)==0)
			cont.push_back(tmp_bits);
#else
		++tmp_cnt;
		tmp_bits=(tmp_bits<<width)+v;
		if(tmp_bits>=bits_thr) {
			cont.push_back(tmp_bits);
			tmp_bits=1;
		}
#endif
	}
	void flush() {
#ifndef ECCL_COMPACT_CODE
		auto rem=tmp_cnt%padded_count;
		if(rem!=0) {
			tmp_cnt+=padded_count-rem;
			cont.push_back(tmp_bits<<(padded_count-rem)*padded_width);
		}
#else
		auto rem=tmp_cnt%nopad_count;
		if(rem!=0) {
			tmp_cnt+=nopad_count-rem;
			cont.push_back(tmp_bits<<(nopad_count-rem)*width);
			tmp_bits=1;
		}
#endif
	};
private:
	Cont& cont;
	std::size_t tmp_cnt=0;
	unsigned int tmp_bits{1};
	static constexpr unsigned int bits_thr=(Type==eccl::seq_type::xna?0x4000'0000:0x4000'0000);
	static constexpr unsigned int width=(Type==eccl::seq_type::xna?3:5);
	static constexpr unsigned int padded_width=(Type==eccl::seq_type::xna?4:8);
	static constexpr unsigned int padded_count=(Type==eccl::seq_type::xna?8:4);
	static constexpr unsigned int nopad_count=(Type==eccl::seq_type::xna?10:6);
};


std::shared_ptr<eccl::fasta_chunk> eccl::fasta_reader::load(eccl::source<char>& str, std::size_t buf_size, unsigned int num_seq) {
	std::shared_ptr<eccl::fasta_chunk> chk{};
	if(st==St::end)
		return chk;
	std::size_t cur_idx{0};
	std::size_t non_seq{0};
	using Coords=fasta_chunk::Coords;
	std::size_t n;
	do {
		while(cur_idx<cur.buf.size()) {
			auto c=cur.buf[cur_idx];
			switch(st) {
			case St::expect_hdr:
				if(c=='>') {
					cur.coords.emplace_back(Coords{cur_idx+1, SIZE_MAX, SIZE_MAX, SIZE_MAX});
					st=St::in_hdr;
					break;
				}
				break;
			case St::in_hdr:
				if((n=cur.buf.find('\n', cur_idx))!=std::string::npos) {
					cur_idx=n;
					auto& coords=cur.coords.back();
					coords.name_e=cur_idx;
					cur.names_total+=coords.name_e-coords.name_s;
					coords.seq_s=cur_idx+1;
					st=St::expect_seq;
					non_seq=0;
					break;
				}
				cur_idx=cur.buf.size()-1;
				break;
			case St::expect_seq:
				if(c=='>') {
					auto& coords=cur.coords.back();
					coords.seq_e=cur_idx;
					cur.seqs_total+=coords.seq_e-coords.seq_s-non_seq;
					st=St::in_hdr;
					if((buf_size>0 && cur.seqs_total>=buf_size) || (num_seq>0 && cur.coords.size()>=num_seq)) {
						Chunk next{};
						next.buf.assign(&cur.buf[cur_idx], &cur.buf[cur.buf.size()]);
						next.coords.emplace_back(Coords{1, SIZE_MAX, SIZE_MAX, SIZE_MAX});
						chk=std::make_shared<Chunk>(std::move(cur));
						cur=std::move(next);
						return chk;
					}
					cur.coords.emplace_back(Coords{cur_idx+1, SIZE_MAX, SIZE_MAX, SIZE_MAX});
					break;
				}
				if(c=='\n') {
					++non_seq;
					break;
				}
				st=St::in_seq;
				break;
			case St::in_seq:
				if((n=cur.buf.find('\n', cur_idx))!=std::string::npos) {
					cur_idx=n;
					++non_seq;
					st=St::expect_seq;
					break;
				}
				cur_idx=cur.buf.size()-1;
				break;
			case St::end:
				assert(0);
			}
			++cur_idx;
		}
		if(hiteof)
			break;
		auto blk=str.get();
		if(blk.eof())
			hiteof=true;
		cur.buf.resize(cur_idx+blk.size());
		std::copy_n(&blk[0], blk.size(), &cur.buf[cur_idx]);
	} while(true);
	switch(st) {
	case St::expect_hdr:
		return {};
	case St::expect_seq:
		break;
	default:
		fprintf(stderr, "wrong st: %d\n", (int)st);
		assert(0);
		throw st;
	}
	st=St::end;
	auto& coords=cur.coords.back();
	coords.seq_e=cur_idx;
	cur.seqs_total+=coords.seq_e-coords.seq_s-non_seq;
	chk=std::make_shared<Chunk>(std::move(cur));
	return chk;
}

template<eccl::seq_type Type>
eccl::seqs<Type>::seqs(const eccl::fasta_chunk& chunk) {
	_b.seq_buf.reserve(chunk.seqs_total/2/sizeof(unsigned int)+chunk.coords.size());
	_b.name_buf.reserve(chunk.names_total);
	_b.coords.reserve(chunk.coords.size());

	detail::BitStream<Type, std::vector<unsigned int>> bits{_b.seq_buf};
	eccl::seq_encoder<Type> b2b;
	_b.max_seq_l=0;
	for(auto coord: chunk.coords) {
		detail::seqs_base::Coords cc;
		while(coord.name_s<coord.name_e)
			_b.name_buf.push_back(chunk.buf[coord.name_s++]);
		cc.name_e=_b.name_buf.size();
		cc.seq_s=_b.seq_buf.size();
		auto sl=coord.seq_e-coord.seq_s;
		while(coord.seq_s<coord.seq_e) {
			auto c=chunk.buf[coord.seq_s++];
			if(c=='\n') {
				--sl;
			} else {
				bits.add(b2b(c).value);
			}
		}
		bits.flush();
		cc.seq_l=sl;
		if(_b.max_seq_l<sl)
			_b.max_seq_l=sl;
		_b.coords.push_back(cc);
	}
}
template eccl::seqs<eccl::seq_type::xna>::seqs(const eccl::fasta_chunk& chunk);
template eccl::seqs<eccl::seq_type::prot>::seqs(const eccl::fasta_chunk& chunk);

eccl::seq_encoder_tables::seq_encoder_tables() noexcept:
_table_xna{0, }, _table_prot{0, }
{
	_rev_xna[0b000]='\x00';
	_rev_xna[0b001]='U';
	_rev_xna[0b010]='C';
	_rev_xna[0b011]='T';
	_rev_xna[0b100]='A';
	_rev_xna[0b101]='G';
	_rev_xna[0b110]='X';
	_rev_xna[0b111]='N';

	_rev_prot[0b00000]='\x00';
	_rev_prot[0b11111]='X';
	_rev_prot[0b11110]='*';
	_rev_prot[0b11101]='-';
	_rev_prot[0b00001]='A';
	_rev_prot[0b00010]='B';
	_rev_prot[0b00011]='C';
	_rev_prot[0b00100]='D';
	_rev_prot[0b00101]='E';
	_rev_prot[0b00110]='F';
	_rev_prot[0b00111]='G';
	_rev_prot[0b01000]='H';
	_rev_prot[0b01001]='I';
	_rev_prot[0b01010]='K';
	_rev_prot[0b01011]='L';
	_rev_prot[0b01100]='M';
	_rev_prot[0b01101]='N';
	_rev_prot[0b01110]='P';
	_rev_prot[0b01111]='Q';
	_rev_prot[0b10000]='R';
	_rev_prot[0b10001]='S';
	_rev_prot[0b10010]='T';
	_rev_prot[0b10011]='U';
	_rev_prot[0b10100]='V';
	_rev_prot[0b10101]='W';
	_rev_prot[0b10110]='Y';
	_rev_prot[0b10111]='Z';

	for(unsigned int i=0; i<_rev_xna.size(); ++i) {
		auto c=_rev_xna[i];
		if(c) {
			_table_xna[std::tolower(c)]=_table_xna[c]=i;
		}
	}
	for(unsigned int i=0; i<_rev_prot.size(); ++i) {
		auto c=_rev_prot[i];
		if(c) {
			_table_prot[std::tolower(c)]=_table_prot[c]=i;
		}
	}
}

template<eccl::seq_type Type>
eccl::seqs_builder<Type>::seqs_builder(std::size_t num_bases, std::size_t names_size, std::size_t num_seqs):
	_b{},
	_str{std::make_shared<detail::BitStream<Type, std::vector<unsigned int>>>(_b.seq_buf)},
	_encoder{}
{
}
template eccl::seqs_builder<eccl::seq_type::xna>::seqs_builder(std::size_t num_bases, std::size_t names_size, std::size_t num_seqs);
template eccl::seqs_builder<eccl::seq_type::prot>::seqs_builder(std::size_t num_bases, std::size_t names_size, std::size_t num_seqs);

template<eccl::seq_type Type, typename T, typename E>
inline static void copy_seq(T& b, E& e, std::string_view seq, int*) {
	for(std::size_t i=0; i<seq.size(); ++i)
		b.add(e(seq[i]).value);
}
template<eccl::seq_type Type, typename T, typename E>
inline static void copy_seq(T& b, E& e, std::string_view seq, char*) {
	static_assert(Type==eccl::seq_type::xna);
	for(std::size_t i=seq.size(); i>0; --i) {
		auto comp=~e(seq[i-1]);
		b.add(comp.value);
	}
}
template<eccl::seq_type Type>
template<bool RC>
void eccl::seqs_builder<Type>::push_back(std::string_view seq) {
	detail::seqs_base::Coords s;
	s.name_e=_b.name_buf.size();
	s.seq_s=_b.seq_buf.size();
	copy_seq<Type>(*_str, _encoder, seq, std::conditional_t<RC, char*, int*>{});
	_str->flush();
	s.seq_l=seq.size();
	if(_b.max_seq_l<seq.size())
		_b.max_seq_l=seq.size();
	_b.coords.push_back(s);
}
template void eccl::seqs_builder<eccl::seq_type::xna>::push_back<true>(std::string_view seq);
template void eccl::seqs_builder<eccl::seq_type::xna>::push_back<false>(std::string_view seq);
template void eccl::seqs_builder<eccl::seq_type::prot>::push_back<false>(std::string_view seq);

#if 0
#endif

eccl::seq_encoder_tables eccl::seq_encoder_tables::_t;

