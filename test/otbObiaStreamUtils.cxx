#include "otbObiaStreamUtils.h"
#include "otbMacro.h"
#include <vector>
#include <array>
#include <algorithm>

int otbObiaStreamUtils(int argc, char * argv[])
{
  std::vector<unsigned long> myVector;
  myVector.push_back(0);
  myVector.push_back(1);
  myVector.push_back(2);
  myVector.push_back(3);

  std::array<float,3> myArray;
  myArray[0]=4.;
  myArray[1]=5.;
  myArray[2]=6.;

  size_t vectorOffset = otb::obia::stream_offset(myVector);

  std::cout<<"Vector offset: "<<vectorOffset<<std::endl;

  double myDouble = 100.;

  size_t arrayOffset = otb::obia::stream_offset(myArray);
  
  std::cout<<"Array offset: "<<arrayOffset<<std::endl;

  size_t doubleOffset = otb::obia::stream_offset(myDouble);

  std::cout<<"Double offset: "<<doubleOffset<<std::endl;

  unsigned short myShort = 12;

  size_t shortOffset = otb::obia::stream_offset(myShort);

  std::cout<<"Unsigned short offset: "<<shortOffset<<std::endl;

  std::vector<char> myStream(vectorOffset+arrayOffset+doubleOffset+shortOffset);
  
  size_t offset = 0;

  otb::obia::to_stream(myStream,myDouble,offset);
  otb::obia::to_stream(myStream,myArray,offset);
  otb::obia::to_stream(myStream,myVector,offset);
  otb::obia::to_stream(myStream,myShort,offset);

  otbControlConditionTestMacro(offset!=myStream.size(),"Offset "<<offset<<" not at end of stream "<<myStream.size());

  offset = 0;
  
  double decodedDouble;
  std::vector<unsigned long> decodedVector;
  std::array<float,3> decodedArray;
  unsigned short decodedShort;

  otb::obia::from_stream(myStream,decodedDouble,offset);
  otb::obia::from_stream(myStream,decodedArray,offset);
  otb::obia::from_stream(myStream,decodedVector,offset);
  otb::obia::from_stream(myStream,decodedShort,offset);

  otbControlConditionTestMacro(offset!=myStream.size(),"Offset not at end of stream.");

  otbControlConditionTestMacro(myDouble != decodedDouble, "Error in double decoding");
  otbControlConditionTestMacro(!std::equal(myArray.begin(), myArray.end(),decodedArray.begin(),decodedArray.end()),"Error in array decoding");
  otbControlConditionTestMacro(!std::equal(myVector.begin(), myVector.end(),decodedVector.begin(),decodedVector.end()),"Error in vector decoding");
  otbControlConditionTestMacro(myShort != decodedShort, "Error in unsigned short decoding");  

  return EXIT_SUCCESS;
}
