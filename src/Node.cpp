#include "Node.h"

Node::Node() :	id(0)
{	
	next = NULL;
}

Node::Node(const std::size_t _id) :	id(_id)
{	
	next = NULL;
}

Node::~Node()
{
	if (next)
		delete next;
}

std::size_t Node::getID() const
{
	return id;
}

Node * Node::getNextNode()
{
	return next;
}

void Node::setID(const std::size_t id)
{
	this->id = id;
}

void Node::setNextNode(Node * next_node)
{
	this->next = next_node;
}
