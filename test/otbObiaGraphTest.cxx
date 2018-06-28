#include "otbObiaGraph.h"

using namespace otb::obia;

int otbObiaGraphTest(int argc, char *argv[])
{

	// Edge tests
	using NodeType = Node< DummyGraphAttribute, DummyGraphAttribute>;
	Edge<NodeType> edge = Edge<NodeType>();
	assert(edge.GetMemorySize() == 8);
	assert(edge.GetNumberOfBytesToSerialize() == 8);
	edge.m_TargetId = 4;
	edge.m_Boundary = 3;

	// Node tests
	NodeType node = NodeType();
	Contour cont1 = Contour();
	cont1.SetStartingCoords(6);
	cont1.AddMove(Contour::Move::RIGHT); cont1.AddMove(Contour::Move::RIGHT); cont1.AddMove(Contour::Move::DOWN); cont1.AddMove(Contour::Move::RIGHT); cont1.AddMove(Contour::Move::DOWN); cont1.AddMove(Contour::Move::LEFT);
	cont1.AddMove(Contour::Move::DOWN); cont1.AddMove(Contour::Move::LEFT); cont1.AddMove(Contour::Move::LEFT); cont1.AddMove(Contour::Move::UP); cont1.AddMove(Contour::Move::UP); cont1.AddMove(Contour::Move::UP);
	node.m_Contour = cont1;

	assert(node.GetFirstPixelCoords() == 6);

	node.SetFirstPixelCoords(5);
	assert(node.GetFirstPixelCoords() == 5);

	auto edge2 = node.AddEdge();
	edge2->m_TargetId = 4;
	edge2->m_Boundary = 3;

	auto lambdaOp = [&](NodeType::EdgeType& edg){
		edg.m_Boundary = 2;
	};
	node.ApplyForEachEdge(lambdaOp);
	assert(node.FindEdge(4)->m_Boundary == 2);

	assert(node.GetMemorySize() == 139);
	assert(node.GetNumberOfBytesToSerialize() == 44);

	// Graph tests

	return EXIT_SUCCESS;
}
