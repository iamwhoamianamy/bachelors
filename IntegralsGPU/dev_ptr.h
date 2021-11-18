#pragma once
#include "cuda_runtime.h"
#include "exeptions.h"

#ifndef DEV_PTR
#define DEV_PTR

namespace cuda_utilities
{

   template <class T> class DevPtr
   {
   private:
      T* _data = nullptr;
      size_t _size = 0;

      DevPtr(const DevPtr&) = delete;
   public:
      DevPtr();
      ~DevPtr<T>();
      DevPtr<T>(size_t size);
      DevPtr<T>(const T* data, size_t size);

      T* Get() const;

      void CopyToHost(T* data);
      void CopyToDevice(const T* data);
      void Init(size_t size);

      size_t Size() const;

      DevPtr& operator=(DevPtr& devPtr);
   };

   template<class T>
   DevPtr<T>::DevPtr()
   {

   }

   template<class T>
   DevPtr<T>::~DevPtr<T>()
   {
      cudaFree(_data);
   }

   template<class T>
   DevPtr<T>::DevPtr<T>(size_t size)
   {
      try
      {
         Init(size);
      }
      catch(Exeption e)
      {
         throw e;
      }
   }

   template<class T>
   DevPtr<T>::DevPtr<T>(const T* data, size_t size)
   {
      try
      {
         Init(size);
         CopyToDevice(data);
      }
      catch(Exeption e)
      {
         throw e;
      }
   }

   template<class T>
   void DevPtr<T>::Init(size_t size)
   {
      cudaFree(_data);
      cudaError_t results = cudaMalloc((void**)&_data,
                                      size * sizeof(T));

      if(results != cudaError_t::cudaSuccess)
         throw MallocExeption();

      _size = size;
   }

   template<class T>
   T* DevPtr<T>::Get() const
   {
      return _data;
   }

   template<class T>
   void DevPtr<T>::CopyToHost(T* data)
   {
      cudaError_t results = cudaMemcpy(data,
                                      _data,
                                      _size * sizeof(T),
                                      cudaMemcpyKind::cudaMemcpyDeviceToHost);

      if(results != cudaError_t::cudaSuccess)
         throw CopyExeption();
   }

   template<class T>
   void DevPtr<T>::CopyToDevice(const T* data)
   {
      cudaError_t results = cudaMemcpy(_data,
                                      data,
                                      _size * sizeof(T),
                                      cudaMemcpyKind::cudaMemcpyHostToDevice);

      if(results != cudaError_t::cudaSuccess)
         throw CopyExeption();
   }

   template<class T>
   size_t DevPtr<T>::Size() const
   {
      return _size;
   }

   template<class T>
   DevPtr<T>& DevPtr<T>::operator=(DevPtr& devPtr)
   {
      cudaFree(_data);

      Init(devPtr._size);
      cudaMemcpy(_data,
                 devPtr._data,
                 devPtr._size * sizeof(T),
                 cudaMemcpyKind::cudaMemcpyDeviceToDevice);

      _size = devPtr._size;

      return *this;
   }
}


#endif // ! DEV_PTR