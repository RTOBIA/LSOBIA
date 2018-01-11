#include "otbObiaStreamUtils.h"

#include <vector>

int otbObiaStreamUtils()
{
  std::vector<unsigned long> myVector;
  myVector.push_back(0);
  myVector.push_back(1);
  myVector.push_back(2);
  myVector.push_back(3);

  size_t vectorOffset = otb::obia::stream_offset(myVector);

  double myDouble = 100.;

  size_t doubleOffset = otb::obia::stream_offset(myDouble);

  unsigned short myShort = 12;

  size_t shortOffset = otb::obia::stream_offset(myShort);

  std::vector<char> myStream(vectorOffset+doubleOffset+shortOffset);
  
  size_t offset = 0;

  otb::obia::to_stream(myStream,myDouble,offset);
  otb::obia::to_stream(myStream,myVector,offset);
  otb::obia::to_stream(myStream,myShort,offset);

  offset = 0;
  
  double decodedDouble;
  std::vector<unsigned long> decodedVector;
  unsigned short decodedShort;

  otb::obia::from_stream(myStream,decodedDouble,offset);
  otb::obia::from_stream(myStream,decodedVector,offset);
  otb::obia::from_stream(myStream,decodedShort,offset);

  return EXIT_SUCCESS;
}
