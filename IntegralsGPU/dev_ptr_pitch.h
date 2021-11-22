#pragma once
#include "cuda_runtime.h"
#include "exeptions.h"

namespace cuda_utilities
{

   template <class T> class DevPtrPitch
   {
   private:
      T* _data = nullptr;
      size_t _size = 0;
      size_t _pitch = 0;
      size_t _width = 0;
      size_t _height = 0;

      DevPtrPitch(const DevPtrPitch&) = delete;
      void InitFields(size_t width, size_t height);

   public:
      DevPtrPitch();
      ~DevPtrPitch<T>();
      DevPtrPitch<T>(size_t width, size_t height);
      DevPtrPitch<T>(const T* data, size_t width, size_t height);

      T* Get() const;

      void CopyToHost(T* data);
      void CopyToDevice(const T* data);
      void Init(size_t width, size_t height);

      size_t Size() const;
      size_t Pitch() const;
      size_t Width() const;
      size_t Height() const;

      DevPtrPitch& operator=(DevPtrPitch& devPtr);
   };

   template<class T>
   inline void DevPtrPitch<T>::InitFields(size_t width, size_t height)
   {
      _width = width;
      _height = height;
      _size = width * height;
   }

   template<class T>
   inline DevPtrPitch<T>::DevPtrPitch()
   {

   }

   template<class T>
   inline DevPtrPitch<T>::~DevPtrPitch<T>()
   {
      cudaFree(_data);
   }

   template<class T>
   inline DevPtrPitch<T>::DevPtrPitch<T>(size_t width, size_t height)
   {
      try
      {
         Init(width, height);

         InitFields(width, height);
      }
      catch(Exeption e)
      {
         throw e;
      }
   }

   template<class T>
   inline DevPtrPitch<T>::DevPtrPitch<T>(const T* data, size_t width, size_t height)
   {
      try
      {
         Init(pitch, width, height);
         CopyToDevice(data);

         InitFields(width, height);
      }
      catch(Exeption e)
      {
         throw e;
      }
   }

   template<class T>
   inline void DevPtrPitch<T>::Init(size_t width, size_t height)
   {
      cudaFree(_data);

      cudaError_t results = cudaMallocPitch((void**)&_data, 
                                            &_pitch, 
                                            width * sizeof(T),
                                            height));

      if(results != cudaError_t::cudaSuccess)
         throw MallocExeption();

      InitFields(width, height);
   }

   template<class T>
   inline T* DevPtrPitch<T>::Get() const
   {
      return _data;
   }

   template<class T>
   inline void DevPtrPitch<T>::CopyToHost(T* data)
   {
      cudaError_t results = cudaMemcpy2D(data, _width * sizeof(T),
                                         _data, _pitch,
                                         _width * sizeof(T), _height,
                                         cudaMemcpyKind::cudaMemcpyDeviceToHost);

      if(results != cudaError_t::cudaSuccess)
         throw CopyExeption();
   }

   template<class T>
   inline void DevPtrPitch<T>::CopyToDevice(const T* data)
   {
      cudaError_t results = cudaMemcpy2D(_data, _pitch,
                                         data, _width * sizeof(T),
                                         _width * sizeof(T), _height,
                                         cudaMemcpyKind::cudaMemcpyHostToDevice);

      if(results != cudaError_t::cudaSuccess)
         throw CopyExeption();
   }

   template<class T>
   inline size_t DevPtrPitch<T>::Size() const
   {
      return _size;
   }

   template<class T>
   inline size_t DevPtrPitch<T>::Pitch() const
   {
      return _pitch;
   }

   template<class T>
   inline size_t DevPtrPitch<T>::Width() const
   {
      return _width;
   }

   template<class T>
   inline size_t DevPtrPitch<T>::Height() const
   {
      return _height;
   }

   template<class T>
   DevPtrPitch<T>& DevPtrPitch<T>::operator=(DevPtrPitch& devPtrPitch)
   {
      cudaFree(_data);

      Init(devPtrPitch._width, devPtrPitch._height);

      cudaMemcpy2D(_data, _pitch,
                   devPtrPitch.data, devPtrPitch._pitch,
                   devPtrPitch._width * sizeof(T), devPtrPitch._height,
                   cudaMemcpyKind::cudaMemcpyDeviceToDevice);

      InitFields(width, height);

      return *this;
   }
}