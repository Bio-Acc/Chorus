#ifndef __ECCL_CUDA_UTILS_HH__
#define __ECCL_CUDA_UTILS_HH__

#include <vector>
#include <cassert>

#include <cuda_runtime.h>

namespace eccl {

enum class device {
	cpu,
	cuda,
};

/*! eg: cudaBlah(...), eccl::check_cuda{"error tag"}; */
class check_cuda {
public:
	explicit constexpr check_cuda(const char* msg) noexcept:
		_msg{msg} { }
private:
	const char* _msg;
	friend void operator,(cudaError_t error, check_cuda checker);
};

void dump_device_info(int device);

//pointer/buffer types:
//
//host_only<t>
//

/*! device raw pointer, no deference in CPU, no memory management */
template<typename T>
class dev_ptr {
public:
	constexpr dev_ptr() noexcept: _p{nullptr} { }
	explicit constexpr dev_ptr(T* p) noexcept: _p{p} { }
	~dev_ptr() { }
	dev_ptr(const dev_ptr&) noexcept =default;
	dev_ptr& operator=(const dev_ptr&) noexcept =default;

	operator T*() const noexcept { return static_cast<T*>(_p); }

private:
	void* _p;
};

template<typename T>
using host_ptr=T*;

namespace detail {

template<typename T>
struct cuda_host_allocator {
	using value_type=T;
	using size_type=std::size_t;
	constexpr cuda_host_allocator() noexcept { }
	~cuda_host_allocator() { }

	T* allocate(size_type n) {
		T* r;
		cudaMallocHost(&r, n*sizeof(T)), check_cuda{"malloc host"};
		//fprintf(stderr, "malloc host: %p %zd\n", r, n*sizeof(T));
		return r;
	}
	void deallocate(T* p, size_type n) {
		//fprintf(stderr, "free.. host: %p %zd\n", p, n*sizeof(T));
		cudaFreeHost(p), check_cuda{"free host"};
	}
};

template<typename T>
struct cuda_dev_allocator {
	using value_type=T;
	using size_type=std::size_t;
	constexpr cuda_dev_allocator() noexcept { }
	~cuda_dev_allocator() { }

	T* allocate(size_type n) {
		T* r;
		cudaMalloc(&r, n*sizeof(T)), check_cuda{"malloc dev"};
		//fprintf(stderr, "malloc dev: %p %zd\n", r, n*sizeof(T));
		//fprintf(stderr, "malloc dev: %p %zd\n", r, n*sizeof(T));
		return r;
	}
	void deallocate(T* p, size_type n) {
		//fprintf(stderr, "free.. dev: %p %zd\n", p, n*sizeof(T));
		cudaFree(p), check_cuda{"free dev"};
	}
};

}

/*! host buffer, accessible to device */
template<typename T>
class host_vector: public std::vector<T, detail::cuda_host_allocator<T>> {
public:
	using std::vector<T, detail::cuda_host_allocator<T>>::vector;
	operator dev_ptr<T>() noexcept {
		return dev_ptr<T>(&(*this)[0]);
	}
};
//using host_vector=

/*! device buffer, std::vector-like API */
template<typename T>
class dev_vector {
public:
	constexpr dev_vector() noexcept: _p{nullptr}, _s{0} { }
	~dev_vector() {
		if(_p) {
			detail::cuda_dev_allocator<T>{}.deallocate(_p, _s);
		}
	}
	dev_vector(const dev_vector&) =delete;
	dev_vector& operator=(const dev_vector&) =delete;
	dev_vector(dev_vector&& r) noexcept: _p{r._p}, _s{r._s} {
		r._p=nullptr;
		r._s=0;
	}
	dev_vector& operator=(dev_vector&& r) noexcept {
		std::swap(_p, r._p);
		std::swap(_s, r._s);
		return *this;
	}

	operator T*() const noexcept { return static_cast<T*>(_p); }
	operator eccl::dev_ptr<T>() const noexcept {
		return eccl::dev_ptr<T>{static_cast<T*>(_p)};
	}
	std::size_t size() const noexcept { return _s; }

	void try_malloc(std::size_t s) {
		if(_p) {
			if(_s<s) {
				free();
				auto ss=_s+std::max<std::size_t>(1024, 16*(s-_s));
				malloc(ss);
			}
		} else {
			malloc(s);
		}
	}
	void malloc(std::size_t s) {
		assert(_p==nullptr);
		_p=detail::cuda_dev_allocator<T>{}.allocate(s);
		_s=s;
	}

	void free() {
		assert(_p);
		detail::cuda_dev_allocator<T>{}.deallocate(_p, _s);
		_p=nullptr;
	}

	template<typename Vec>
			void try_malloc_h2d(const Vec& vec, cudaStream_t stream) {
				try_malloc(vec.size());
				h2d(vec, stream);
			}
	void malloc_h2d(const host_vector<T>& vec, cudaStream_t stream) {
		malloc(vec.size());
		h2d(vec, stream);
	}
	template<typename Alloc>
	void h2d(const std::vector<T, Alloc>& vec, cudaStream_t stream) {
		h2d(&vec[0], vec.size(), stream);
	}
	void h2d(const T* p, std::size_t s, cudaStream_t stream) {
		assert(s<=_s);
		//fprintf(stderr, "H %p -> D %p %zd\n", &_q[0], _p, _q.size()*sizeof(T));
		cudaMemcpyAsync(_p, p, s*sizeof(T), cudaMemcpyHostToDevice, stream), check_cuda{"memcpy async to device"};
	}
	void h2d(const T* p, std::size_t s) {
		assert(s<=_s);
		cudaMemcpy(_p, p, s*sizeof(T), cudaMemcpyHostToDevice), check_cuda{"memcpy to device"};
	}
	void d2h(T* p, std::size_t s) {
		assert(s<=_s);
		//fprintf(stderr, "D2h %p %p %zd\n", _p, p, s*sizeof(T));
		cudaMemcpy(p, _p, s*sizeof(T), cudaMemcpyDeviceToHost), check_cuda{"memcpy to host"};
	}
	void d2h(T* p, std::size_t s, cudaStream_t stream) {
		assert(s<=_s);
		//fprintf(stderr, "D2h async %p %p %zd\n", _p, p, s*sizeof(T));
		cudaMemcpyAsync(p, _p, s*sizeof(T), cudaMemcpyDeviceToHost, stream), check_cuda{"memcpy async to host"};
	}
	template<typename Alloc>
	void d2h(std::vector<T, Alloc>& vec, cudaStream_t stream) {
		d2h(&vec[0], vec.size(), stream);
	}

private:
	T* _p;
	std::size_t _s{0};
};

/////////////////////////////////////////

class event {
public:
	constexpr event() noexcept: _e{nullptr} { }
	~event() {
		if(_e)
			cudaEventDestroy(_e), check_cuda{"create event"};
	}
	event(const event&) =delete;
	event& operator=(const event&) =delete;
	//event(event&& r) =delete;
	//event& operator=(event&& r) =delete;

	void create() {
		cudaEventCreate(&_e), check_cuda{"create event"};
	}
	void record(cudaStream_t stream) {
		cudaEventRecord(_e, stream), check_cuda{"record event"};
	}

	operator cudaEvent_t() const noexcept { return _e; }

private:
	cudaEvent_t _e;
};


class stream {
public:
	constexpr stream() noexcept: _s{nullptr} { }
	~stream() {
		if(_s)
			cudaStreamDestroy(_s), check_cuda{"destroy stream"};
	}
	stream(const stream&) =delete;
	stream& operator=(const stream&) =delete;
	//stream(stream&& r) =delete;
	//stream& operator=(stream&& r) =delete;

	void create() {
		cudaStreamCreate(&_s), check_cuda{"create stream"};
	}

	operator cudaStream_t() const noexcept { return _s; }

	void wait(cudaEvent_t event) {
		cudaStreamWaitEvent(_s, event, cudaEventWaitDefault), check_cuda{"wait event"};
	}
	void sync() {
		cudaStreamSynchronize(_s), check_cuda{"sync stream"};
	}

private:
	cudaStream_t _s;
};

}

#endif
