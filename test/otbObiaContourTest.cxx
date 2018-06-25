#include "otbObiaContour.h"

using namespace otb::obia;

int otbObiaContourTest(int argc, char *argv[])
{
	// Test constructor
	Contour cont1 = Contour();
	Contour cont2 = Contour();

	// Test StartingCoords accessors
	cont1.SetStartingCoords(6);
	cont2.SetStartingCoords(3);
	assert(cont1.GetStartingCoords() == 6);

	// Test AddMoves
	cont1.AddMove(Contour::Move::RIGHT); cont1.AddMove(Contour::Move::RIGHT); cont1.AddMove(Contour::Move::DOWN); cont1.AddMove(Contour::Move::RIGHT); cont1.AddMove(Contour::Move::DOWN); cont1.AddMove(Contour::Move::LEFT);
	cont1.AddMove(Contour::Move::DOWN); cont1.AddMove(Contour::Move::LEFT); cont1.AddMove(Contour::Move::LEFT); cont1.AddMove(Contour::Move::UP); cont1.AddMove(Contour::Move::UP); cont1.AddMove(Contour::Move::UP);
	cont2.AddMove(Contour::Move::RIGHT); cont2.AddMove(Contour::Move::DOWN); cont2.AddMove(Contour::Move::RIGHT); cont2.AddMove(Contour::Move::DOWN); cont2.AddMove(Contour::Move::LEFT); cont2.AddMove(Contour::Move::DOWN);
	cont2.AddMove(Contour::Move::LEFT); cont2.AddMove(Contour::Move::UP); cont2.AddMove(Contour::Move::LEFT); cont2.AddMove(Contour::Move::UP); cont2.AddMove(Contour::Move::RIGHT); cont2.AddMove(Contour::Move::UP);

	// Test Iterator
	Contour::ContourIterator it = cont1.Begin();
	assert(it.GetMove() == Contour::Move::RIGHT);
	it.GoToNext(); it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::DOWN);
	assert(it.IsAtEnd() == false);
	it.GoToNext(); it.GoToNext(); it.GoToNext(); it.GoToNext(); it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::UP);
	assert(it.IsAtEnd() == true);

	// Test GetMemorySize
	assert(cont1.GetMemorySize() == sizeof(Contour)+12*UInt8Size);

	// Test serialization
	uint64_t dummy = 0;
	std::vector<char> serialCont1;
	serialCont1.assign(cont1.GetNumberOfBytesToSerialize(), char());
	cont1.Serialize(serialCont1, dummy);
	Contour cont1b = Contour();
	dummy = 0;
	cont1b.DeSerialize(serialCont1, dummy);
	assert(cont1b.GetStartingCoords() == 6);
	it = cont1b.Begin();
	assert(it.GetMove() == Contour::Move::RIGHT);
	it.GoToNext(); it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::DOWN);

	// Test MergeWith
	cont1.MergeWith(cont2, 5, 4);
	assert(cont1.GetStartingCoords() == 3);
	it = cont1.Begin();
	assert(it.GetMove() == Contour::Move::RIGHT);
	it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::DOWN);
	it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::LEFT);
	it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::LEFT);
	it.GoToNext(); it.GoToNext(); it.GoToNext();
	assert(it.GetMove() == Contour::Move::UP);

	// Test GenerateBorderPixels
	Contour::CoordsSet borderCells = Contour::CoordsSet();
	cont2.GenerateBorderPixels(borderCells, 5);
	auto it2 = borderCells.find(3);
	assert(it2 != borderCells.end());
	it2 = borderCells.find(7);
	assert(it2 != borderCells.end());
	it2 = borderCells.find(8);
	assert(it2 != borderCells.end());
	it2 = borderCells.find(9);
	assert(it2 != borderCells.end());
	it2 = borderCells.find(13);
	assert(it2 != borderCells.end());


	return EXIT_SUCCESS;
}

