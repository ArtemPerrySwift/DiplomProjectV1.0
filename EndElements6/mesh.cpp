#include "mesh.h"

bool operator >> (ifstream& in, MeshData& mData) { return mData.readData(in); }
