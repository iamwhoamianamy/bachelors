#pragma once
#include "cuda_runtime.h"
#include "cuda_helper.hpp"
#include "math.hpp"

namespace cuda
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

      T* data() const;

      void copyToHost(T* result);
      void copyToDevice(const T* source);
      void allocateSpaceForData(size_t size, size_t padding);

      size_t size() const;
      size_t sizePadded() const;

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
         allocateSpaceForData(size, padding);
      }
      catch(std::exception e)
      {
         throw e;
      }
   }

   template<class T>
   DevPtr<T>::DevPtr<T>(const T* data, size_t size, size_t padding)
   {
      try
      {
         allocateSpaceForData(size, padding);
         copyToDevice(data);
      }
      catch(std::exception e)
      {
         throw e;
      }
   }

   template<class T>
   void DevPtr<T>::allocateSpaceForData(size_t size, size_t padding)
   {
      cudaFree(_data);

      _sizePadded = math::nextDevisible(size, padding);
      _size = size;

      cudaError_t result = cudaMalloc((void**)&_data,
                                      _sizePadded * sizeof(T));

      if(result != cudaError_t::cudaSuccess)
         throw std::exception("deb");
   }

   template<class T>
   T* DevPtr<T>::data() const
   {
      return _data;
   }

   template<class T>
   void DevPtr<T>::copyToHost(T* result)
   {
      cudaError_t resultMessage = 
         cudaMemcpy(result, _data, _size * sizeof(T),
                    cudaMemcpyKind::cudaMemcpyDeviceToHost);

      if(resultMessage != cudaError_t::cudaSuccess)
         throw std::exception();
   }

   template<class T>
   void DevPtr<T>::copyToDevice(const T* source)
   {
      cudaError_t results = cudaMemcpy(_data,
                                      source,
                                      _size * sizeof(T),
                                      cudaMemcpyKind::cudaMemcpyHostToDevice);

      if(results != cudaError_t::cudaSuccess)
         throw std::exception();
   }

   template<class T>
   size_t DevPtr<T>::size() const
   {
      return _size;
   }

   template<class T>
   inline size_t DevPtr<T>::sizePadded() const
   {
      return _sizePadded;
   }

   template<class T>
   DevPtr<T>& DevPtr<T>::operator=(DevPtr& devPtr)
   {
      allocateSpaceForData(devPtr._size);

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