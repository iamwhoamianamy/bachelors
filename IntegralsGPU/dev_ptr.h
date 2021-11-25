#pragma once
#include "cuda_runtime.h"
#include "cuda_helper.h"
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
      size_t _sizePadded = 0;

      DevPtr(const DevPtr&) = delete;
   public:
      DevPtr();
      ~DevPtr<T>();
      DevPtr<T>(size_t size, size_t padding = 0);
      DevPtr<T>(const T* data, size_t size, size_t padding = 0);

      T* Get() const;

      void CopyToHost(T* data);
      void CopyToDevice(const T* data);
      void Init(size_t size, size_t padding);

      size_t Size() const;
      size_t SizePadded() const;

      DevPtr& operator=(DevPtr& devPtr);
      DevPtr& operator=(DevPtr&& devPtr);
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
   DevPtr<T>::DevPtr<T>(size_t size, size_t padding)
   {
      try
      {
         Init(size, padding);
      }
      catch(Exeption e)
      {
         throw e;
      }
   }

   template<class T>
   DevPtr<T>::DevPtr<T>(const T* data, size_t size, size_t padding)
   {
      try
      {
         Init(size, padding);
         CopyToDevice(data);
      }
      catch(Exeption e)
      {
         throw e;
      }
   }

   template<class T>
   void DevPtr<T>::Init(size_t size, size_t padding)
   {
      cudaFree(_data);

      _sizePadded = nextDevisible(size, padding);
      _size = size;

      cudaError_t result = cudaMalloc((void**)&_data,
                                      _sizePadded * sizeof(T));

      if(result != cudaError_t::cudaSuccess)
         throw MallocExeption();
   }

   template<class T>
   T* DevPtr<T>::Get() const
   {
      return _data;
   }

   template<class T>
   void DevPtr<T>::CopyToHost(T* data)
   {
      cudaError_t result = cudaMemcpy(data,
                                      _data,
                                      _size * sizeof(T),
                                      cudaMemcpyKind::cudaMemcpyDeviceToHost);

      if(result != cudaError_t::cudaSuccess)
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
   inline size_t DevPtr<T>::SizePadded() const
   {
      return _sizePadded;
   }

   template<class T>
   DevPtr<T>& DevPtr<T>::operator=(DevPtr& devPtr)
   {
      Init(devPtr._size);

      cudaMemcpy(_data,
                 devPtr._data,
                 devPtr._sizePadded * sizeof(T),
                 cudaMemcpyKind::cudaMemcpyDeviceToDevice);

      _size = devPtr._size;
      _sizePadded = devPtr._sizePadded;

      return *this;
   }

   template<class T>
   DevPtr<T>& DevPtr<T>::operator=(DevPtr&& devPtr)
   {
      _data = devPtr._data;
      _size = devPtr._size;
      _sizePadded = devPtr._sizePadded;

      devPtr._data = nullptr;
      devPtr._size = 0;
      devPtr._sizePadded = 0;

      return *this;
   }
}


#endif // ! DEV_PTR