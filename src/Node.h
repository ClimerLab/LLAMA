#ifndef NODE_H
#define NODE_H

#include <cstddef>

class Node
{
public:
	Node();
	Node(const std::size_t);
	~Node();
	std::size_t getID() const;
	Node *getNextNode();
	void setID(const std::size_t);
	void setNextNode(Node *);

private:
	std::size_t id;
	Node *next;
};

#endif