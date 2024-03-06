#ifndef __ECCL_CORE_HH__
#define __ECCL_CORE_HH__

#include <memory>
#include <string>
#include <future>
#include <cstring>
#include <ostream>
#include <cassert>

/*! core types/definitions/functions, cheap to include */


#ifdef __CUDACC__
#define CUDA_HOST __host__
#define CUDA_DEVICE __device__
#else
#define CUDA_HOST
#define CUDA_DEVICE
#endif

/*! 10 necleotides per unsigned int (instead of 8) */
#define ECCL_COMPACT_CODE

namespace eccl {

enum class seq_type {
	xna, /*! DNA or RNA */
	//dna,
	//rna,
	prot,
};

template<seq_type Type>
struct code {
	unsigned char value{0};
	constexpr static unsigned int width=(Type==seq_type::xna?3:5);
#ifndef ECCL_COMPACT_CODE
	constexpr static unsigned int n_per_word=(Type==seq_type::xna?8:4);
#else
	constexpr static unsigned int n_per_word=(Type==seq_type::xna?10:6);
#endif
	/*! complement neucleotide */
	code<Type> operator~() const noexcept {
		static_assert(Type==seq_type::xna);
		return code<Type>{static_cast<unsigned char>((~value)&0b111)};
	}
};
using nucleotide=code<seq_type::xna>;
using amino_acid=code<seq_type::prot>;

template<seq_type Type>
struct pair {
	unsigned int value{0};
	CUDA_HOST CUDA_DEVICE constexpr pair(code<Type> a, code<Type> b) noexcept:
		value{(static_cast<unsigned int>(a.value)<<a.width)+b.value} { }
};

template<seq_type Type>
inline unsigned int padded_len(unsigned int len) noexcept {
#ifndef ECCL_COMPACT_CODE
	return Type==seq_type::xna?(len+7)/8*8:(len+3)/4*4;
#else
	//return (len+9)/10*10;
	return Type==seq_type::xna?(len+39)/40*40:(len+23)/24*24;
#endif
}

template<eccl::seq_type Type>
inline CUDA_DEVICE code<Type> get_code(const unsigned int* buf, unsigned long idx);
template<>
inline CUDA_DEVICE code<eccl::seq_type::xna> get_code(const unsigned int* buf, unsigned long idx) {
#ifndef ECCL_COMPACT_CODE
	auto v=buf[idx/8];
	return {static_cast<unsigned char>((v>>((7-idx%8)*4))&0x0f)};
#else
	auto v=buf[idx/10];
	return {static_cast<unsigned char>((v>>((9-idx%10)*3))&0b0111)};
#endif
};
template<>
inline CUDA_DEVICE code<eccl::seq_type::prot> get_code(const unsigned int* buf, unsigned long idx) {
#ifndef ECCL_COMPACT_CODE
	auto v=buf[idx/4];
	return {static_cast<unsigned char>((v>>((3-idx%4)*8))&0x1f)};
#else
	auto v=buf[idx/6];
	return {static_cast<unsigned char>((v>>((5-idx%6)*5))&0b011111)};
#endif
};

template<typename T>
class chunk {
public:
	constexpr chunk() noexcept: _p{nullptr}, _s{0} { }
	explicit chunk(std::size_t size): _p{std::make_unique<T[]>(size)}, _s{size<<1} { }

	explicit operator bool() const noexcept { return _p.get(); }
	std::size_t size() const noexcept { return _s>>1; }
	T& operator[](std::size_t i) noexcept { return _p[i]; }
	const T& operator[](std::size_t i) const noexcept { return _p[i]; }
	bool eof() const noexcept { return _s&1; }

	void shrink(std::size_t size, bool eof=false) noexcept {
		assert(size<=(_s>>1));
		_s=(size<<1)|(eof?1:0);
	}

private:
	std::unique_ptr<T[]> _p;
	/*! use last bit for eof */
	std::size_t _s;
};

template<typename T>
class source {
public:
	chunk<T> get() {
		auto chk=_fut.get();
		_fut=prepare();
		return chk;
	}

protected:
	constexpr source() noexcept { }
	void post_ctor() {
		_fut=prepare();
	}
	virtual ~source() { }

	virtual std::future<chunk<T>> prepare() =0;
private:
	std::future<chunk<T>> _fut;
};


}

#if __cplusplus < 201703L
namespace std {
class string_view {
public:
	constexpr string_view(const char* p, std::size_t s) noexcept:
		_p{p}, _s{s} { }
	string_view(const char* p):
		_p{p}, _s{std::strlen(p)} { }
	std::size_t size() const noexcept { return _s; }
	const char& operator[](std::size_t i) const noexcept { return _p[i]; }
private:
	const char* _p;
	std::size_t _s;
};
inline bool operator!=(const std::string& a, std::string_view b) noexcept {
	return a.compare(0, a.size(), &b[0], b.size())!=0;
}
inline bool operator==(const std::string& a, std::string_view b) noexcept {
	return !(a!=b);
}
inline std::ostream& operator<<(std::ostream& oss, std::string_view b) {
	return oss.write(&b[0], b.size());
}
}
#endif

#endif
