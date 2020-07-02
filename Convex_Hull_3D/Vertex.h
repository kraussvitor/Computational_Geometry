#ifndef VERTEX_H_INCLUDED
#define VERTEX_H_INCLUDED

#include <Eigen/Dense>

class Vertex
{
	private:
		Eigen::Vector3d location;

	public:
		Vertex();
		Vertex(Eigen::Vector3d location);
		Eigen::Vector3d getLocation() const;

};

#endif // !
