#include "dev_ptr.h"
#include "cuda_runtime.h"

//template<class T>
//DevPtr<T>::DevPtr()
//{
//
//}

//template<class T>
//DevPtr<T>::DevPtr(const DevPtr& devPtr)
//{
//   cudaFree(_data);
//   _data = devPtr.Get();
//   _size = devPtr.Size();
//}


//
//template<class T>
//DevPtr<T>& DevPtr<T>::operator=(DevPtr<T>& devPtr)
//{
//   cudaFree(_data);
//   _data = devPtr.Get();
//   _size = devPtr.Size();
//
//   return this;
//}
//
//template<class T>
//T& DevPtr<T>::operator[](int i)
//{
//   return _data[i];
//}
//
//template<class T>
//const T& DevPtr<T>::operator[](int i) const
//{
//   return _data[i];
//}
//
//template<class T>
//void DevPtr<T>::CopyFromHost(const T* data)
//{
//   cudaError_t result = cudaMemcpy(static_cast<void*>(_data),
//                                   static_cast<void*>(data),
//                                   _size * sizeof(T),
//                                   cudaMemcpyKind::cudaMemcpyHostToDevice);
//   if(result != cudaError_t::cudaSuccess)
//      throw CopyError();
//}
//
//template<class T>
//void DevPtr<T>::CopyToHost(T* data)
//{
//   cudaError_t result = cudaMemcpy((void**)_data,
//                                   (void**)data,
//                                   _size * sizeof(T),
//                                   cudaMemcpyKind::cudaMemcpyHostToDevice);
//   if(result != cudaError_t::cudaSuccess)
//      throw CopyError();
//}




