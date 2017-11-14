#ifndef __otbObiaStreamUtils_h
#define __otbObiaStreamUtils_h

#include <vector>
#include <array>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <typeinfo>

namespace otb
{
namespace obia
{


  template <typename T> uint64_t stream_offset(const T& value)
    {
      return sizeof(T);
    }

  template <typename T> uint64_t stream_offset(const std::vector<T>& value)
    {
      return sizeof(value.size()) + value.size() * sizeof(T);
    }

  template <typename T> void to_stream(std::vector<char> & stream, const T& value, uint64_t& position)
    {
      uint64_t offset = stream_offset(value);
      assert(position+offset<=stream.size());
      std::memcpy(&(stream[position]),&value,offset);
      //std::cout<<"Wrote "<<typeid(T).name()<<" from "<<position<<" to "<<position+offset<<std::endl;
      position+=offset;      
    }

  template <typename T> void to_stream(std::vector<char> & stream, const std::vector<T>& value, uint64_t& position)
    {
      uint64_t offset = stream_offset(value);
      auto nb_elements = value.size();

      /* if(position+offset>stream.size()) */
      /* 	std::cerr<<sizeof(T)<<" "<<nb_elements<<" "<<position<<" "<<offset<<" "<<stream.size()<<std::endl; */

      assert(position+offset<=stream.size());
      // First copy vector size
      std::memcpy(&(stream[position]),&nb_elements,sizeof(nb_elements));

      if(nb_elements>0)
	{
	  // Then copy vector values
	  std::memcpy(&(stream[position+sizeof(nb_elements)]),&(value[0]),offset-sizeof(nb_elements));
	}

      //std::cout<<"Wrote a vector of "<<nb_elements<<" of "<<typeid(T).name()<<" from "<<position<<" to "<<position+offset<<std::endl;
      position+=offset;      
    }  

  template<typename T> void from_stream(const std::vector<char> & stream, T& value, uint64_t& position)
    {
      uint64_t offset = stream_offset(value);

      assert(position+offset<=stream.size());

      std::memcpy(&value,&(stream[position]),offset);
      //std::cout<<"Read a "<<typeid(T).name()<<" from "<<position<<" to "<<position+offset<<std::endl;
      position+=offset;
    }

 
  template<typename T> void from_stream(const std::vector<char> & stream, std::vector<T>& value, uint64_t& position)
    {
      
      auto nb_elements = value.size();

      // Decode nb elements first
      std::memcpy(&nb_elements,&(stream[position]),sizeof(nb_elements));

      // Resize container
      value.resize(nb_elements);
      uint64_t offset = stream_offset(value);

      // Compute offset
      assert(position+offset<=stream.size());

      if(nb_elements > 0)
	{
	  // Decode vector elements
	  std::memcpy(&(value[0]),&(stream[position+sizeof(nb_elements)]),offset-sizeof(nb_elements));
        }
      //std::cout<<"Read a vector of "<<nb_elements<<" of "<<typeid(T).name()<<" from "<<position<<" to "<<position+offset<<std::endl;
      position+=offset;
    }

 template<typename T> void from_stream(const std::vector<char> & stream, T& value)
    {
      uint64_t dummy = 0;
      from_stream(stream,value,dummy);
    }

}
}

#endif
