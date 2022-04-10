#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <list>
#include <type_traits>
#include "proglogact.h"

using namespace std;

struct MeshData
{
	virtual bool readData(ifstream& in) = 0;
};

bool operator >> (ifstream& in, MeshData& mData); 

struct Border: MeshData
{
	int fun_num;
	int vertL, vertR;
	int horizL, horizR;
	int heightL, heightR;
	
	Border(ifstream& in) { readData(in);}

	Border() {};

	virtual bool readData(ifstream &in)
	{
		in >> fun_num >> horizL >> horizR >> vertL >> vertR >> heightL >> heightR;
		return true;
	}

};

struct BorderFirstCond : Border
{
	BorderFirstCond(ifstream& in) : Border(in){}
	BorderFirstCond(){};
	//Структура создана для возможного расширения
};

struct BorderSecondCond : Border
{
	BorderSecondCond(ifstream& in) : Border(in) {}
	BorderSecondCond(){};
	//Структура создана для возможного расширения
};

struct BorderThirdCond : Border
{
	double betta;
	BorderThirdCond(ifstream& in) : Border(in) { in >> betta; }
	BorderThirdCond() {};

	bool readData(ifstream &in) override
	{
		Border::readData(in);
		in >> betta;
		return true;
	}
};

//template <class StorageTypeCond, class = enable_if_t<is_cond_type<StorageTypeCond>>>

template<class TypeMeshData>
constexpr bool is_cond_type = is_base_of<MeshData, TypeMeshData>::value;

template <class TypeMeshData>
//template <class StorageTypeCond, class = enable_if<is_cond_type<StorageTypeCond>>>
struct StoreMeshData
{
	static_assert(is_cond_type<TypeMeshData>, L"Задан неверный тип для хранилища данных сетки");
	TypeMeshData* data;
	int count;

	//StoreMeshData() {};
	StoreMeshData(int count = 0) 
	{
		if (count < 0) { writeErr("Указанное количество границ < 0 "); return; }
		this->count = count;

		data = new TypeMeshData[count];
	};

	virtual bool readStorage(ifstream& in, string errCase = "Ошибка считывания из файла ", int count = 0)
	{
		if (count == 0)
		{
			try { in >> count; }
			catch (...) { writeErr(errCase); return false; }
		}
			
		if(count < 0) { writeErr("Указанное количество границ < 0 "); return false; }

		this->count = count;

		data = new TypeMeshData[count];

		try
		{
			for (int i = 0; i < count; i++)
				data[i].readData(in);
				//in >> data[i];
		}
		catch (...)
		{
			writeErr(errCase);
			return false;
		}
		return true;
	}
};

struct Area: MeshData // параметры подобласти
{
	double lambda;
	double gamma;
	int lXW, rXW; // Номера левой и правой границы области по X
	int lYW, rYW; // Номера левой и правой границы области по Y
	int lZW, rZW; // Номера левой и правой границы области по Z
	bool readData(ifstream& in) override
	{
		try 
		{
			in >> lambda >> gamma;
			in >> lXW >> rXW;
			in >> lYW >> rYW;
			in >> lZW >> rZW; 
		}
		catch (...) { writeErr("Ошибка при чтении параметров областей"); return false; }

		return true;
	}
};

struct SepParam : MeshData // параметры разбиения отрезка сетки
{
	int n; // количество отрезков разбиения
	double q; // коэффициент разбиения
	bool readData(ifstream& in) override
	{
		try { in >> n >> q; }
		catch (...) { writeErr("Ошибка при чтении параметров разбиения"); return false; }
		return true;
	}
};

struct Coord : MeshData
{
	double x;
	double y;
	double z;

	Coord(){x = y = z = 0;}
	Coord(double x, double y, double z) 
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void write(ofstream& out)
	{
		out << x << " " << y << " " << z;
	}

	bool readData(ifstream& in) override
	{
		try { in >> x >> y >> z; }
		catch (...) { writeErr("Ошибка при чтении координат"); return false; }
		return true;
	}

	Coord operator +(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x + coord2.x;
		newCoord.y = this->y + coord2.y;
		newCoord.z = this->z + coord2.z;
		return newCoord;
	}

	Coord operator -(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x - coord2.x;
		newCoord.y = this->y - coord2.y;
		newCoord.z = this->z - coord2.z;
		return newCoord;
	}

	Coord operator /(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x / coord2.x;
		newCoord.y = this->y / coord2.y;
		newCoord.z = this->z / coord2.z;
		return newCoord;
	}

	Coord operator *(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x * coord2.x;
		newCoord.y = this->y * coord2.y;
		newCoord.z = this->z * coord2.z;
		return newCoord;
	}

	Coord operator /(const double& a)
	{
		Coord newCoord;
		newCoord.x = this->x / a;
		newCoord.y = this->x / a;
		newCoord.z = this->x / a;
		return newCoord;
	}

	Coord operator *(const double& a)
	{
		Coord newCoord;
		newCoord.x = this->x * a;
		newCoord.y = this->x * a;
		newCoord.z = this->x * a;
		return newCoord;
	}

	Coord& operator +=(const Coord& coord2)
	{
		Coord newCoord;
		this->x += coord2.x;
		this->y += coord2.y;
		this->z += coord2.z;
		return *this;
	}

	Coord& operator *=(const Coord& coord2)
	{
		Coord newCoord;
		this->x *= coord2.x;
		this->y *= coord2.y;
		this->z *= coord2.z;
		return *this;
	}

	Coord& operator *=(const double& a)
	{
		Coord newCoord;
		this->x *= a;
		this->y *= a;
		this->z *= a;
		return *this;
	}

};

struct StoreMeshKnots : StoreMeshData<Coord>
{
	int nX, nY, nZ;
	int nKnots;
	int knotsInLine, knotsInPlane;

	StoreMeshKnots(int nX, int nY, int nZ){initStore(nX, nY, nZ);}

	Coord& operator [](int i) { return data[i]; }

	bool initStore(int nX, int nY, int nZ)
	{
		this->nX = nX;
		this->nY = nY;
		this->nZ = nZ;
		knotsInLine = nX;
		knotsInPlane = nY * knotsInLine;
		nKnots = nX * nY * nZ;
		StoreMeshData<Coord>::StoreMeshData<Coord>(nKnots);
		return true;
	}

	bool readStorage(ifstream& in, string errCase = "Ошибка считывания сетки из файла ", int n = 0) override
	{
		in >> nX >> nY >> nZ;
		knotsInLine = nX;
		knotsInPlane = nY * knotsInLine;
		nKnots = nX * nY * nZ;
		nKnots = nX * nY * nZ;
		StoreMeshData<Coord>::readStorage(in, errCase, nKnots);
		return true;
	}

	void writeAllX(char seperator = ' ')
	{
		for (int i = 0; i < nKnots; i++)
			cout << data[i].x << " ";
		cout << endl << endl;
	}

	void writeAllY(char seperator = ' ')
	{
		for (int i = 0; i < nKnots; i++)
			cout << data[i].y << " ";
		cout << endl << endl;
	}

	void writeAllZ(char seperator = ' ')
	{
		for (int i = 0; i < nKnots; i++)
			cout << data[i].z << " ";
		cout << endl << endl;
	}

	int getKnotIndex(int ix, int jy, int kz)
	{
		return kz* knotsInPlane + jy * knotsInLine + ix;
	}

	int writeMesh(string fileName)
	{
		int i, j, k, p;
		ofstream out;
		out.open(fileName);
		out << nX << " " << nY << " " << nZ << endl;
		for (k = 0, p = 0; k < nZ; k++)
		{
			for (j = 0; j < nY; j++)
			{
				for (i = 0; i < nX; i++, p++)
				{
					data[p].write(out);
					out << " ";
				}
					
				out << endl;
			}
			out << endl;
		}
		return 0;
	}
	StoreMeshKnots() {};
};

struct ContainerSepParam
{
	StoreMeshData<SepParam> sepX;
	StoreMeshData<SepParam> sepY;
	StoreMeshData<SepParam> sepZ;

	int nSepX, nSepY, nSepZ;

	ContainerSepParam(int nSepX, int nSepY, int nSepZ, ifstream& in) { readSepParams(nSepX, nSepY, nSepZ, in); }

	ContainerSepParam() {};

	bool readSepParams(int nSepX, int nSepY, int nSepZ, ifstream& in)
	{
		this->nSepX = nSepX;
		sepX.readStorage(in, "Ошибка чтения параметров разбиения по оси x из файла", nSepX);

		this->nSepY = nSepY;
		sepY.readStorage(in, "Ошибка чтения параметров разбиения по оси y из файла", nSepY);

		this->nSepZ = nSepZ;
		sepZ.readStorage(in, "Ошибка чтения параметров разбиения по оси z из файла", nSepZ);
		
		return true;
	}
};

struct ContainerIXYZW
{
	int* IXW, * IYW, * IZW;
	int nIXW, nIYW, nIZW;

	ContainerIXYZW(ContainerSepParam containerSepParam) { create(containerSepParam); }
	ContainerIXYZW() {};

	void create(ContainerSepParam containerSepParam)
	{
		countICoordW(containerSepParam.sepX, IXW, nIXW);
		countICoordW(containerSepParam.sepY, IYW, nIYW);
		countICoordW(containerSepParam.sepZ, IZW, nIZW);;
	}

	void writeIX(char seperator = ' ')
	{
		for (int i = 0; i < nIXW; i++)
			cout << IXW[i] << seperator;
		cout << endl << endl;
	}

	void writeIY(char seperator = ' ')
	{
		for (int i = 0; i < nIYW; i++)
			cout << IYW[i] << seperator;
		cout << endl << endl;
	}

	void writeIZ(char seperator = ' ')
	{
		for (int i = 0; i < nIZW; i++)
			cout << IZW[i] << seperator;
		cout << endl << endl;
	}

private:
	void countICoordW(StoreMeshData<SepParam> sepCoord, int*& ICoordW, int& nICoordW)
	{
		SepParam* sepParamC = sepCoord.data;
		
		int nSepC = sepCoord.count;
		ICoordW = new int[nSepC + 1];
		int nC = 0;
		for (int i = 0; i < nSepC; i++)
		{
			ICoordW[i] = nC;
			nC += sepParamC[i].n;
		}
		ICoordW[nSepC] = nC;
		nICoordW = nSepC + 1;
	}
};

struct StoreCrushedMeshKnots : StoreMeshKnots
{
	ContainerIXYZW IXYZW;

	StoreCrushedMeshKnots() {};
	StoreCrushedMeshKnots(StoreMeshKnots XYZW, ContainerSepParam containerSepParam) { initDividesStore(XYZW, containerSepParam); }

	void initDividesStore(StoreMeshKnots XYZW, ContainerSepParam containerSepParam)
	{
		IXYZW.create(containerSepParam);
		this->initStore(IXYZW.IXW[IXYZW.nIXW - 1] + 1, IXYZW.IYW[IXYZW.nIYW - 1] + 1, IXYZW.IZW[IXYZW.nIZW - 1] + 1);

		knotsInLine = nX; // Количество точек в горизонтальной линии
		knotsInPlane = nY * knotsInLine; // Количество точек в вертикальной линии
		nKnots = nZ * knotsInPlane;

		int knotsWInLine = XYZW.nX; // Количество основных точек в горизонтальной линии
		int knotsWInPlane = XYZW.nY * knotsWInLine; // Количество основных точек в вертикальной линии

		SepParam* sepX = containerSepParam.sepX.data;
		SepParam* sepY = containerSepParam.sepY.data;
		SepParam* sepZ = containerSepParam.sepZ.data;

		SepParam sepParamX, sepParamY, sepParamZ;

		/*
		double qHoriz, qVert, qHeight;
		double nHoriz, nVert, nHeight;

		double koef;

		double dx, dy, dz;
		double x, y, z;
		int i, p;
		
		int pWCurr, // текущая точка (основные координаты)
			pWNextHeight, // следующая точка по высоте (основные координаты)
			pWNextHoriz, // следующая горизонтальная точка (основные координаты)
			pWNextVert; // следующая вертикальная точка (основные координаты)
		*/

		int pCurr, // текущая точка (общие координаты)
			pNextHeight, // следующая точка по высоте (общие координаты)
			pNextHoriz, // следующая горизонтальная точка (общие координаты)
			pNextVert; // следующая вертикальная точка (общие координаты)

		
		Coord dCoord;
		Coord coord;
		
		int iX = 0, iY = 0, iZ = 0;

		int nXW = XYZW.nX,
			nYW = XYZW.nY,
			nZW = XYZW.nZ;

		int* IXW = IXYZW.IXW;
		int* IYW = IXYZW.IYW;
		int* IZW = IXYZW.IZW;

		int iXW, iYW, iZW;

		copyBasePoints(XYZW);

		for (iZW = 0; iZW < nZW; iZW++)
		{
			iZ = IZW[iZW];
			/*sepParamZ = sepZ[iZW];
			
			nHeight = sepZ[iZW].n;
			qHeight = sepZ[iZW].q;
			iY = 0;
			*/
			for (iYW = 0; iYW < nYW; iYW++)
			{
				iY = IYW[iYW];
				/*sepParamY = sepY[iYW];	
				/*
				nVert = sepY[iYW].n;
				qVert = sepY[iYW].q;
				iX = 0;
				*/
				for (iXW = 0; iXW < nXW - 1; iXW++)
				{
					sepParamX = sepX[iXW];
					iX = IXW[iXW];

					pCurr = getKnotIndex(iX, iY, iZ);
					//pNextHoriz = getKnotIndex(iX + 1, iY, iZ);
					pNextHoriz = pCurr + sepX[iXW].n;

					copyLineOfBasePoints(sepParamX, XYZW, pNextHoriz, pCurr);
					/*
					pWCurr       = XYZW.getKnotIndex(iXW + 0, iYW + 0, iZW + 0);
					pWNextHoriz  = XYZW.getKnotIndex(iXW + 1, iYW + 0, iZW + 0);
					pWNextVert   = XYZW.getKnotIndex(iXW + 0, iYW + 1, iZW + 0); 
					pWNextHeight = XYZW.getKnotIndex(iXW + 0, iYW + 0, iZW + 1);
					
					/*
					pWCurr = iZW * knotsWInPlane + iYW * knotsWInLine + iXW;
					pWNextHoriz = pWCurr + 1;
					pWNextVert = pWCurr + knotsWInLine;
					pWNextHeight = pWCurr + knotsWInPlane;
					
					
					pCurr		= this->getKnotIndex(iX + 0, iY + 0, iZ + 0);
					pNextHoriz  = this->getKnotIndex(iX + 1, iY + 0, iZ + 0);
					pNextVert   = this->getKnotIndex(iX + 0, iY + 1, iZ + 0);
					pNextHeight = this->getKnotIndex(iX + 0, iY + 0, iZ + 1);

					/*
					pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
					pNextHoriz = pCurr + sepX[iXW].n;
					pNextVert = pCurr + sepY[iYW].n * knotsInLine;
					pNextHeight = pCurr + sepZ[iZW].n * knotsInPlane;
					
					
					/*
					nHoriz = sepX[iXW].n;
					qHoriz = sepX[iXW].q;
					

					/*Проходим в направлении горизонтальной ломанной
					if (iXW != nXW - 1)
					{
						copyLineOfBasePoints(sepParamX, pWNextHoriz, pWCurr, pNextHoriz, pCurr, XYZW);
						/*
						if (qHoriz == 1)
						{
							dCoord = (XYZW[pWNextHoriz] - XYZW[pWCurr]) / nHoriz;
							/*
							dx = (XW[pWNextHoriz] - XW[pWCurr]) / nHoriz;
							dy = (YW[pWNextHoriz] - YW[pWCurr]) / nHoriz;
							dz = (ZW[pWNextHoriz] - ZW[pWCurr]) / nHoriz;
							
						}
						else
						{
							koef = (1 - qHoriz) / (1 - pow(qHoriz, nHoriz));
							dCoord = (XYZW[pWNextHoriz] - XYZW[pWCurr]) * koef;
							/*
							dx = (XW[pWNextHoriz] - XW[pWCurr]) * koef;
							dy = (YW[pWNextHoriz] - YW[pWCurr]) * koef;
							dz = (ZW[pWNextHoriz] - ZW[pWCurr]) * koef;
							
						}

						coord = XYZW[pWCurr];
						/*
						x = XW[pWCurr];
						y = YW[pWCurr];
						z = ZW[pWCurr];
						
						for (p = pCurr; p < pNextHoriz; p++)
						{
							XYZ[p] = coord;
							/*
							X[p] = x;
							Y[p] = y;
							Z[p] = z;
							

							coord += dCoord;
							/*
							x += dx;
							y += dy;
							z += dz;
							

							dCoord *= qHoriz;
							/*
							dx *= qHoriz;
							dy *= qHoriz;
							dz *= qHoriz;
							
						}
						
					}

					Проходим в направлении вертикальной ломанной
					if (iYW != nYW - 1)
					{
						copyLineOfBasePoints(sepParamY, pWNextVert, pWCurr, pNextVert, pCurr, XYZW);
						/*
						if (qVert == 1)
						{
							dCoord = (XYZW[pWNextHoriz] - XYZW[pWCurr]) / nVert;
							/*
							dx = (XW[pWNextVert] - XW[pWCurr]) / nVert;
							dy = (YW[pWNextVert] - YW[pWCurr]) / nVert;
							dz = (ZW[pWNextVert] - ZW[pWCurr]) / nVert;
							
						}
						else
						{
							koef = (1 - qVert) / (1 - pow(qVert, nVert));
							dCoord = (XYZW[pWNextHoriz] - XYZW[pWCurr]) * koef;
							/*
							dx = (XW[pWNextVert] - XW[pWCurr]) * koef;
							dy = (YW[pWNextVert] - YW[pWCurr]) * koef;
							dz = (ZW[pWNextVert] - ZW[pWCurr]) * koef;
							
						}

						coord = XYZW[pWCurr];
						/*
						x = XW[pWCurr];
						y = YW[pWCurr];
						z = ZW[pWCurr];
						

						for (p = pCurr; p < pNextVert; p += knotsInLine)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qVert;
							dy *= qVert;
							dz *= qVert;
						}
						
					}


					Проходим в направлении ломанной по высоте
					if (iZW != nZW - 1)
					{
						copyLineOfBasePoints(sepParamZ, pWNextHeight, pWCurr, pNextHeight, pCurr, XYZW);
						/*
						if (qHeight == 1)
						{
							dx = (XW[pWNextHeight] - XW[pWCurr]) / nHeight;
							dy = (YW[pWNextHeight] - YW[pWCurr]) / nHeight;
							dz = (ZW[pWNextHeight] - ZW[pWCurr]) / nHeight;
						}
						else
						{
							koef = (1 - qHeight) / (1 - pow(qHeight, nHeight));
							dx = (XW[pWNextHeight] - XW[pWCurr]) * koef;
							dy = (YW[pWNextHeight] - YW[pWCurr]) * koef;
							dz = (ZW[pWNextHeight] - ZW[pWCurr]) * koef;
						}

						x = XW[pWCurr];
						y = YW[pWCurr];
						z = ZW[pWCurr];

						for (p = pCurr; p < pNextHeight; p += knotsInPlane)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qHeight;
							dy *= qHeight;
							dz *= qHeight;
						}
						
					}
					*/
				}
			}
		}
		/*XYZ[nKnots - 1] = XYZW[XYZW.nKnots - 1];
		
		X[nKnots - 1] = XW[nW - 1];
		Y[nKnots - 1] = YW[nW - 1];
		Z[nKnots - 1] = ZW[nW - 1];
		*/

		for (iZW = 0; iZW < nZW; iZW++)
		{
			iZ = IZW[iZW];
			for (iYW = 0; iYW < nYW - 1; iYW++)
			{
				iY = IYW[iYW];
				sepParamY = sepY[iYW];
				for (iXW = 0; iXW < nXW - 1; iXW++)
				{
					for (iX = IXW[iXW] + 1; iX < IXW[iXW + 1] + 1; iX++)
					{
						/*
						nHoriz = sepX[iXW].n;
						qHoriz = sepX[iXW].q;
						iX = IXW[iXW];
						*/
						pCurr = getKnotIndex(iX, iY, iZ);
						//pNextVert = getKnotIndex(iX, iY + 1, iZ);
						pNextVert = pCurr + sepParamY.n * knotsInLine;

						copyLineOfBasePoints(sepParamX, XYZW, pNextVert, pCurr, knotsInLine);
						/*
						pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
						pNextVert = pCurr + sepParamY.n * knotsInLine;

						if (qHoriz == 1)
						{
							dx = (X[pNextHoriz] - X[pCurr]) / nHoriz;
							dy = (Y[pNextHoriz] - Y[pCurr]) / nHoriz;
							dz = (Z[pNextHoriz] - Z[pCurr]) / nHoriz;
						}
						else
						{
							koef = (1 - qHoriz) / (1 - pow(qHoriz, nHoriz));
							dx = (X[pNextHoriz] - X[pCurr]) * koef;
							dy = (Y[pNextHoriz] - Y[pCurr]) * koef;
							dz = (Z[pNextHoriz] - Z[pCurr]) * koef;
						}

						x = X[pCurr] + dx;
						y = Y[pCurr] + dy;
						z = Z[pCurr] + dz;

						dx *= qHoriz;
						dy *= qHoriz;
						dz *= qHoriz;

						for (p = pCurr + 1; p < pNextHoriz; p++)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qHoriz;
							dy *= qHoriz;
							dz *= qHoriz;
						}
						*/
					}
				}
			}
		}
		//int j;

		for (iZW = 0; iZW < nZW - 1; iZW++)
		{
			iZ = IZW[iZW];
			sepParamZ = sepZ[iZW];
			/*
			nHeight = sepZ[iZW].n;
			qHeight = sepZ[iZW].q;
			//j = 0;
			//pCurr = iZ * knotsInPlane;
			//pNextHeight = pCurr + nHeight * knotsInPlane;
			*/
			for (iYW = 0; iYW < nYW - 1; iYW++)
			{
				for (iY = IYW[iYW]; iY < IYW[iYW + 1] + 1; iY++)
				{
					for (iXW = 0; iXW < nXW - 1; iXW++)
					{
						for (iX = IXW[iXW]; iX < IXW[iXW + 1] + 1; iX++)
						{
							if (iX == IXW[iXW] && iY == IYW[iYW]) continue;

							pCurr = getKnotIndex(iX, iY, iZ);
							//pNextHeight = getKnotIndex(iX, iY, iZ + 1);
							pNextHeight = pCurr + sepParamZ.n * knotsInPlane;

							copyLineOfBasePoints(sepParamZ, XYZW, pNextHeight, pCurr, knotsInPlane);
							/*
							pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
							pNextHeight = pCurr + nHeight * knotsInPlane;
							if (qHeight == 1)
							{
								dx = (X[pNextHeight] - X[pCurr]) / nHeight;
								dy = (Y[pNextHeight] - Y[pCurr]) / nHeight;
								dz = (Z[pNextHeight] - Z[pCurr]) / nHeight;
							}
							else
							{
								koef = (1 - qHeight) / (1 - pow(qHeight, nHeight));
								dx = (X[pNextHeight] - X[pCurr]) * koef;
								dy = (Y[pNextHeight] - Y[pCurr]) * koef;
								dz = (Z[pNextHeight] - Z[pCurr]) * koef;
							}

							x = X[pCurr] + dx;
							y = Y[pCurr] + dy;
							z = Z[pCurr] + dz;

							dx *= qHeight;
							dy *= qHeight;
							dz *= qHeight;

							for (p = pCurr + knotsInPlane; p < pNextHeight; p += knotsInPlane)
							{
								X[p] = x;
								Y[p] = y;
								Z[p] = z;

								x += dx;
								y += dy;
								z += dz;

								dx *= qHeight;
								dy *= qHeight;
								dz *= qHeight;
							}
							*/
						}
					}
				}
			}
		}
	}

private: 
	void copyBasePoints(StoreMeshKnots& XYZW)
	{
		Coord* XYZ = this->data;
		int iXW, iYW, iZW;
		int* IXW = IXYZW.IXW;
		int* IYW = IXYZW.IYW;
		int* IZW = IXYZW.IZW;
		int nXW = XYZW.nX,
			nYW = XYZW.nY,
			nZW = XYZW.nZ;

		int pCurr, pWCurr;
		for (iZW = 0; iZW < nXW; iZW++)
		{
			for (iYW = 0; iYW < nYW; iYW++)
			{
				for (iXW = 0; iXW < nXW; iXW++)
				{
					pWCurr = XYZW.getKnotIndex(iXW, iYW, iZW);
					pCurr = getKnotIndex(IXW[iXW], IYW[iYW], IZW[iZW]);

					XYZ[pCurr] = XYZW[pWCurr];
				}
			}
		}
	}

	void copyLineOfBasePoints(SepParam sepParam, StoreMeshKnots& XYZW, int pNext, int pCurr, int pStep = 1)
	{
		Coord* XYZ = this->data;
		Coord coord, dCoord;
		double koef;
		double q = sepParam.q;
		int n = sepParam.n;
		if (q == 1)
		{
			dCoord = (XYZ[pNext] - XYZ[pCurr]) / n;
			/*
			dx = (XW[pWNextHoriz] - XW[pWCurr]) / nHoriz;
			dy = (YW[pWNextHoriz] - YW[pWCurr]) / nHoriz;
			dz = (ZW[pWNextHoriz] - ZW[pWCurr]) / nHoriz;
			*/
		}
		else
		{
			koef = (1 - q) / (1 - pow(q, n));
			dCoord = (XYZ[pNext] - XYZ[pCurr]) * koef;
				/*
				dx = (XW[pWNextHoriz] - XW[pWCurr]) * koef;
				dy = (YW[pWNextHoriz] - YW[pWCurr]) * koef;
				dz = (ZW[pWNextHoriz] - ZW[pWCurr]) * koef;
				*/
		}

		/*
		x = XW[pWCurr];
		y = YW[pWCurr];
		z = ZW[pWCurr];
		*/

		coord = XYZ[pCurr] + dCoord;		
		dCoord *= q;

		for (int p = pCurr + pStep; p < pNext; p += pStep)
		{
			XYZ[p] = coord;
				/*
				X[p] = x;
				Y[p] = y;
				Z[p] = z;
				*/

			coord += dCoord;
				/*
				x += dx;
				y += dy;
				z += dz;
				*/

			dCoord *= q;
				/*
				dx *= qHoriz;
				dy *= qHoriz;
				dz *= qHoriz;
				*/
		}
	}
};

struct ContainerBorders
{
	/*
	BorderFirstCond* firstCondAr;
	BorderSecondCond* secondCondAr;
	BorderThirdCond* thirdCondAr;
	*/
	StoreMeshData<BorderFirstCond> firstCondStor;
	StoreMeshData<BorderSecondCond> secondCondStor;
	StoreMeshData<BorderThirdCond> thirdCondStor;
	//BorderStorCondType<int> hey;
	/*
	int nFirstCond;
	int nSecondCond;
	int nThirdCond;
	*/
	bool readBordersCond1(ifstream& in) { return firstCondStor.readStorage(in); }

	bool readBordersCond2(ifstream& in) { return secondCondStor.readStorage(in); }

	bool readBordersCond3(ifstream& in) { return thirdCondStor.readStorage(in); }

	bool readBorders(ifstream& in) { return  readBordersCond1(in) && readBordersCond2(in) && readBordersCond3(in); }

	bool readBorders(string fileName) 
	{ 
		ifstream in;
		in.open(fileName);
		bool res = readBordersCond1(in) && readBordersCond2(in) && readBordersCond3(in);
		in.close();
		return res; 
	}
};

struct CalculationArea // пространственная сетка
{
	StoreMeshKnots XYZW;
	StoreCrushedMeshKnots XYZ;
	StoreMeshData<Area> Areas;
	ContainerBorders borders;
	/*
	double* XW, *YW, *ZW; // координаты опорных узлов
	int nW; // общее количество опорных узлов
	int nXW, nYW, nZW; // количество опорных ломанных по соответствующим осям

	int* IXW, *IYW, *IZW;

	double *X, *Y, *Z; // координаты всех узлов сетки
	int n; // общее количество узлов
	int nX, nY, nZ; // общее количество ломанных по соответствущим осям
	
	int nArea; // количество подобластей
	Area* areas; // подобласти
	*/
	
	CalculationArea(string fileNameCord, string fileNameSep, string fileNameBorders, bool writeInCordFile = false)
	{
		ifstream in(fileNameCord);
		
		/*Считывание основных ломанных*/
		XYZW.readStorage(in);
		/*
		in >> nXW >> nYW >> nZW; // считывание количества основынх ломанных по сотвествующим осям
		
		nW = nXW * nYW * nZW; 

		XW = new double[nW];
		YW = new double[nW];
		ZW = new double[nW];

		int i;
		for (i = 0; i < nW; i++)
		{
			in >> XW[i] >> YW[i] >> ZW[i];
		}
		*/

		/*Считывание информации о подобластях*/
		Areas.readStorage(in);
		/*
		in >> nArea;
		areas = new Area[nArea];
		for (i = 0; i < nArea; i++)
		{
			in >> areas[i].lambda >> areas[i].gamma;
			in >> areas[i].lXW >> areas[i].rXW;
			in >> areas[i].lYW >> areas[i].rYW;
			in >> areas[i].lZW >> areas[i].rZW;
		}
		in.close();
		*/

		cout << "XW" << endl;
		XYZW.writeAllX();

		cout << "YW" << endl;
		XYZW.writeAllY();

		cout << "ZW" << endl;
		XYZW.writeAllZ();
		/*
		for (i = 0; i < nW; i++)
			cout << XW[i] << " ";
		cout << endl << endl;

		
		for (i = 0; i < nW; i++)
			cout << YW[i] << " ";
		cout << endl << endl;

		
		for (i = 0; i < nW; i++)
			cout << ZW[i] << " ";
		cout << endl << endl;
		*/
		in.close();
		/*Считывание параметров разбиения*/
		in.open(fileNameSep);

		ContainerSepParam containerSepParam;
		containerSepParam.readSepParams(XYZW.nX - 1, XYZW.nY - 1, XYZW.nZ - 1, in);

		/*
		ContainerIXYZW IXYZW(containerSepParam);
		XYZ.initStore(IXYZW.IXW[IXYZW.nIXW - 1] + 1, IXYZW.IYW[IXYZW.nIYW - 1] + 1, IXYZW.IZW[IXYZW.nIZW - 1] + 1);
		
		int nSepX = nXW - 1, 
			nSepY = nYW - 1, 
			nSepZ = nZW - 1;
		
		SepParam* sepX = new SepParam[nSepX],
				* sepY = new SepParam[nSepY],
				* sepZ = new SepParam[nSepZ];

		
		IXW = new int[nXW];
		IYW = new int[nYW];
		IZW = new int[nZW];

		nX = 0;
		for (i = 0; i < nSepX; i++)
		{
			IXW[i] = nX;
			in >> sepX[i].n >> sepX[i].q;
			nX += sepX[i].n;
		}
		IXW[i] = nX;
		nX += 1;

		cout << "IXW" << endl;
		for (i = 0; i < nXW; i++)
			cout << IXW[i] << " ";
		cout << endl << endl;


		nY = 0;
		for (i = 0; i < nSepY; i++)
		{
			IYW[i] = nY;
			in >> sepY[i].n >> sepY[i].q;
			nY += sepY[i].n;
		}
		IYW[i] = nY;
		nY += 1;

		cout << "IYW" << endl;
		for (i = 0; i < nYW; i++)
			cout << IYW[i] << " ";
		cout << endl << endl;

		nZ = 0;
		for (i = 0; i < nSepZ; i++)
		{
			IZW[i] = nZ;
			in >> sepZ[i].n >> sepZ[i].q;
			nZ += sepZ[i].n;
		}
		IZW[i] = nZ;
		nZ += 1;

		cout << "IZW" << endl;
		for (i = 0; i < nZW; i++)
			cout << IZW[i] << " ";
		cout << endl << endl;
		

		int knotsInLine = XYZ.nX; // Количество точек в горизонтальной линии
		int knotsInPlane = XYZ.nY * knotsInLine; // Количество точек в вертикальной линии
		int knotsN = XYZ.nZ * knotsInPlane;

		int knotsWInLine = XYZW.nX; // Количество основных точек в горизонтальной линии
		int knotsWInPlane = XYZW.nY * knotsWInLine; // Количество основных точек в вертикальной линии
		*/
		XYZ.initDividesStore(XYZW, containerSepParam);
		//StoreCrushedMeshKnots(XYZW, containerSepParam);

		borders.readBorders(fileNameBorders);
		/*
		double heigtHorizPx, heigtHorizPy, heigtHorizPz;
		double heigtVertPx, heigtVertPy, heigtVertPz;
		double heigtPx, heigtPy, heigtPz;

		double heigtdxHorizPx, heigtdyHorizPy, heigtdzHorizPz;
		double heigtdxVertPx, heigdytVertPy, heigtdzVertPz;
		double heigtdxPx, heigtdyPy, heigtdzPz;

		double heigtHorizPdx, heigtHorizPdy, heigtHorizPdz;
		double heigtVertPdx, heigtVertPdy, heigtVertPdz;
		double heigtPdx, heigtPdy, heigtPdz;
		*/
		in.close();
		/*Генерация сетки*/	
		//int iXW, iYW, iZW;
		//*
		//int inX, inY, inZ;
		//double iqX, iqY, iqZ;
		//int ix, iy, iz;
		//*/
		//
		//X = new double[nKnots];
		//Y = new double[nKnots];
		//Z = new double[nKnots];
		//*/
		//
		//int iKnot = 0;
		//int iKnotPrev = 0;
		//int c, l;
		//
		//*
		//double dxHoriz, dyHoriz, dzHoriz;
		//double dxHeght, dyHeight, dzHeight;
		//double dxVert, dyVert, dzVert;
		//*/
		//double qHoriz, qVert, qHeight;
		//double nHoriz, nVert, nHeight;
		//
		//double koef;
		//
		//int pWCurr, // текущая точка (основные координаты)
		//	pWNextHeight, // следующая точка по высоте (основные координаты)
		//	pWNextHoriz, // следующая горизонтальная точка (основные координаты)
		//	pWNextVert; // следующая вертикальная точка (основные координаты)
		//
		//int pCurr, // текущая точка (общие координаты)
		//	pNextHeight, // следующая точка по высоте (общие координаты)
		//	pNextHoriz, // следующая горизонтальная точка (общие координаты)
		//	pNextVert; // следующая вертикальная точка (общие координаты)
		//
		//double dx, dy, dz;
		//double x, y, z;
		//int iX = 0, iY = 0, iZ = 0;
		//int p;
		//int i;
		//
		//int nXW = XYZW.nX, 
		//	nYW = XYZW.nY,
		//	nZW = XYZW.nZ;
		//
		//for (i = 0, iZW = 0; iZW < nZW; iZW++)
		//{
		//	nHeight = sepZ[iZW].n;
		//	qHeight = sepZ[iZW].q;
		//	iZ = IZW[iZW];
		//	iY = 0;
		//	for (iYW = 0; iYW < nYW; iYW++)
		//	{
		//		nVert = sepY[iYW].n;
		//		qVert = sepY[iYW].q;
		//		iY = IYW[iYW];
		//		iX = 0;
		//		for (iXW = 0; iXW < nXW; iXW++)
		//		{
		//			pWCurr = iZW * knotsWInPlane + iYW * knotsWInLine + iXW;
		//			pWNextHoriz = pWCurr + 1;
		//			pWNextVert = pWCurr + knotsWInLine;
		//			pWNextHeight = pWCurr + knotsWInPlane;
		//
		//			iX = IXW[iXW];
		//			pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
		//			pNextHoriz = pCurr + sepX[iXW].n;
		//			pNextVert = pCurr + sepY[iYW].n * knotsInLine;
		//			pNextHeight = pCurr + sepZ[iZW].n * knotsInPlane;
		//
		//			/*Проходим в направлении горизонтальной ломанной*/
		//			if (iXW != nXW - 1)
		//			{
		//				nHoriz = sepX[iXW].n;
		//				qHoriz = sepX[iXW].q;
		//				if (qHoriz == 1)
		//				{
		//					dx = (XW[pWNextHoriz] - XW[pWCurr]) / nHoriz;
		//					dy = (YW[pWNextHoriz] - YW[pWCurr]) / nHoriz;
		//					dz = (ZW[pWNextHoriz] - ZW[pWCurr]) / nHoriz;
		//				}
		//				else
		//				{
		//					koef = (1 - qHoriz) / (1 - pow(qHoriz, nHoriz));
		//					dx = (XW[pWNextHoriz] - XW[pWCurr]) * koef;
		//					dy = (YW[pWNextHoriz] - YW[pWCurr]) * koef;
		//					dz = (ZW[pWNextHoriz] - ZW[pWCurr]) * koef;
		//				}
		//
		//				x = XW[pWCurr];
		//				y = YW[pWCurr];
		//				z = ZW[pWCurr];
		//
		//				for (p = pCurr; p < pNextHoriz; p++)
		//				{
		//					X[p] = x;
		//					Y[p] = y;
		//					Z[p] = z;
		//
		//					x += dx;
		//					y += dy;
		//					z += dz;
		//
		//					dx *= qHoriz;
		//					dy *= qHoriz;
		//					dz *= qHoriz;
		//				}
		//			}
		//
		//			/*Проходим в направлении вертикальной ломанной*/
		//			if (iYW != nYW - 1)
		//			{
		//				if (qVert == 1)
		//				{
		//					dx = (XW[pWNextVert] - XW[pWCurr]) / nVert;
		//					dy = (YW[pWNextVert] - YW[pWCurr]) / nVert;
		//					dz = (ZW[pWNextVert] - ZW[pWCurr]) / nVert;
		//				}
		//				else
		//				{
		//					koef = (1 - qVert) / (1 - pow(qVert, nVert));
		//					dx = (XW[pWNextVert] - XW[pWCurr]) * koef;
		//					dy = (YW[pWNextVert] - YW[pWCurr]) * koef;
		//					dz = (ZW[pWNextVert] - ZW[pWCurr]) * koef;
		//				}
		//
		//				x = XW[pWCurr];
		//				y = YW[pWCurr];
		//				z = ZW[pWCurr];
		//
		//				for (p = pCurr; p < pNextVert; p += knotsInLine)
		//				{
		//					X[p] = x;
		//					Y[p] = y;
		//					Z[p] = z;
		//
		//					x += dx;
		//					y += dy;
		//					z += dz;
		//
		//					dx *= qVert;
		//					dy *= qVert;
		//					dz *= qVert;
		//				}
		//			}
		//			
		//
		//			/*Проходим в направлении ломанной по высоте*/
		//			if (iZW != nZW - 1)
		//			{
		//				if (qHeight == 1)
		//				{
		//					dx = (XW[pWNextHeight] - XW[pWCurr]) / nHeight;
		//					dy = (YW[pWNextHeight] - YW[pWCurr]) / nHeight;
		//					dz = (ZW[pWNextHeight] - ZW[pWCurr]) / nHeight;
		//				}
		//				else
		//				{
		//					koef = (1 - qHeight) / (1 - pow(qHeight, nHeight));
		//					dx = (XW[pWNextHeight] - XW[pWCurr]) * koef;
		//					dy = (YW[pWNextHeight] - YW[pWCurr]) * koef;
		//					dz = (ZW[pWNextHeight] - ZW[pWCurr]) * koef;
		//				}
		//
		//				x = XW[pWCurr];
		//				y = YW[pWCurr];
		//				z = ZW[pWCurr];
		//
		//				for (p = pCurr; p < pNextHeight; p += knotsInPlane)
		//				{
		//					X[p] = x;
		//					Y[p] = y;
		//					Z[p] = z;
		//
		//					x += dx;
		//					y += dy;
		//					z += dz;
		//
		//					dx *= qHeight;
		//					dy *= qHeight;
		//					dz *= qHeight;
		//				}
		//			}
		//			
		//		}
		//	}
		//}
		//X[knotsN - 1] = XW[nW - 1];
		//Y[knotsN - 1] = YW[nW - 1];
		//Z[knotsN - 1] = ZW[nW - 1];
		//
		//for (i = 0, iZW = 0; iZW < nZW; iZW++)
		//{
		//	iZ = IZW[iZW];
		//	for (iYW = 0; iYW < nYW - 1; iYW++)
		//	{
		//		for (iY = IYW[iYW] + 1; iY < IYW[iYW + 1]; iY++)
		//		{
		//			for (iXW = 0; iXW < nXW - 1; iXW++)
		//			{
		//				nHoriz = sepX[iXW].n;
		//				qHoriz = sepX[iXW].q;
		//				iX = IXW[iXW];
		//				pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
		//				pNextHoriz = pCurr + sepX[iXW].n;
		//
		//				if (qHoriz == 1)
		//				{
		//					dx = (X[pNextHoriz] - X[pCurr]) / nHoriz;
		//					dy = (Y[pNextHoriz] - Y[pCurr]) / nHoriz;
		//					dz = (Z[pNextHoriz] - Z[pCurr]) / nHoriz;
		//				}
		//				else
		//				{
		//					koef = (1 - qHoriz) / (1 - pow(qHoriz, nHoriz));
		//					dx = (X[pNextHoriz] - X[pCurr]) * koef;
		//					dy = (Y[pNextHoriz] - Y[pCurr]) * koef;
		//					dz = (Z[pNextHoriz] - Z[pCurr]) * koef;
		//				}
		//
		//				x = X[pCurr] + dx;
		//				y = Y[pCurr] + dy;
		//				z = Z[pCurr] + dz;
		//
		//				dx *= qHoriz;
		//				dy *= qHoriz;
		//				dz *= qHoriz;
		//
		//				for (p = pCurr + 1; p < pNextHoriz; p++)
		//				{
		//					X[p] = x;
		//					Y[p] = y;
		//					Z[p] = z;
		//
		//					x += dx;
		//					y += dy;
		//					z += dz;
		//
		//					dx *= qHoriz;
		//					dy *= qHoriz;
		//					dz *= qHoriz;
		//				}
		//
		//			}
		//		}
		//	}
		//}
		//int j;
		//
		//for (i = 0, iZW = 0; iZW < nZW - 1; iZW++)
		//{
		//	iZ = IZW[iZW];
		//	nHeight = sepZ[iZW].n;
		//	qHeight = sepZ[iZW].q;
		//	j = 0;
		//	pCurr = iZ * knotsInPlane;
		//	pNextHeight = pCurr + nHeight*knotsInPlane;
		//	
		//	for (iYW = 0; iYW < nYW - 1; iYW++)
		//	{
		//		for (iY = IYW[iYW]; iY < IYW[iYW + 1] + 1; iY++)
		//		{
		//			for (iXW = 0; iXW < nXW - 1; iXW++)
		//			{
		//				for (iX = IXW[iXW]; iX < IXW[iXW + 1] + 1; iX++, pNextHeight++)
		//				{
		//					if (iX == IXW[iXW] && iY == IYW[iYW]) continue;
		//					pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
		//					pNextHeight = pCurr + nHeight*knotsInPlane;
		//					if (qHeight == 1)
		//					{
		//						dx = (X[pNextHeight] - X[pCurr]) / nHeight;
		//						dy = (Y[pNextHeight] - Y[pCurr]) / nHeight;
		//						dz = (Z[pNextHeight] - Z[pCurr]) / nHeight;
		//					}
		//					else
		//					{
		//						koef = (1 - qHeight) / (1 - pow(qHeight, nHeight));
		//						dx = (X[pNextHeight] - X[pCurr]) * koef;
		//						dy = (Y[pNextHeight] - Y[pCurr]) * koef;
		//						dz = (Z[pNextHeight] - Z[pCurr]) * koef;
		//					}
		//
		//					x = X[pCurr] + dx;
		//					y = Y[pCurr] + dy;
		//					z = Z[pCurr] + dz;
		//
		//					dx *= qHeight;
		//					dy *= qHeight;
		//					dz *= qHeight;
		//
		//					for (p = pCurr + knotsInPlane; p < pNextHeight; p += knotsInPlane)
		//					{
		//						X[p] = x;
		//						Y[p] = y;
		//						Z[p] = z;
		//
		//						x += dx;
		//						y += dy;
		//						z += dz;
		//
		//						dx *= qHeight;
		//						dy *= qHeight;
		//						dz *= qHeight;
		//					}
		//
		//				}
		//			}
		//		}
		//	}
		//}
	}
	/*
	int readBorders(string fileName)
	{
		ifstream in;
		in.open(fileName);

		//ContainerBorders borderStorage;
		borders.readBorders(in);
		/*
		BorderFirstCond borderFirstCond_buf;
		BorderSecondCond  borderSecondCond_buf;
		BorderThirdCond borderThirdCond_buf;

		BorderFirstCond* borderFirstCondAr;
		BorderSecondCond* borderSecondCondAr;
		BorderThirdCond* borderThirdCondAr;

		int nFirstCond;
		in >> nFirstCond;
		borders.nFirstCond = nFirstCond;
		borders.firstCondAr = new BorderFirstCond[nFirstCond];
		borderFirstCondAr = borders.firstCondAr;
		//borders = new border[N_bord];
		for (int i = 0; i < nFirstCond; i++)
		{
			in >> borderFirstCond_buf.fun_num;
			in >> borderFirstCond_buf.horizL;
			in >> borderFirstCond_buf.horizR;
			in >> borderFirstCond_buf.vertL;
			in >> borderFirstCond_buf.vertR;
			in >> borderFirstCond_buf.heightL;
			in >> borderFirstCond_buf.heightR;
			borderFirstCondAr[i] = borderFirstCond_buf;
		}

		int nSecondCond;
		in >> nSecondCond;
		borders.nSecondCond = nSecondCond;
		if (nSecondCond == 0)
			borders.secondCondAr = NULL;
		else
			borders.secondCondAr = new BorderSecondCond[nFirstCond];
		borderSecondCondAr = borders.secondCondAr;
		//borders = new border[N_bord];
		for (int i = 0; i < nSecondCond; i++)
		{
			//in >> borderSecondCond_buf.fun_num;
			in >> borderSecondCond_buf.fun_num;
			in >> borderSecondCond_buf.horizL;
			in >> borderSecondCond_buf.horizR;
			in >> borderSecondCond_buf.vertL;
			in >> borderSecondCond_buf.vertR;
			in >> borderSecondCond_buf.heightL;
			in >> borderSecondCond_buf.heightR;
			borderSecondCondAr[i] = borderSecondCond_buf;
		}

		int nThirdCond;
		in >> nThirdCond;
		borders.nThirdCond = nThirdCond;
		if (nThirdCond == 0)
			borders.thirdCondAr = NULL;
		else
			borders.thirdCondAr = new BorderThirdCond[nFirstCond];
		borderThirdCondAr = borders.thirdCondAr;
		//borders = new border[N_bord];
		for (int i = 0; i < nThirdCond; i++)
		{
			in >> borderThirdCond_buf.fun_num;
			in >> borderThirdCond_buf.horizL;
			in >> borderThirdCond_buf.horizR;
			in >> borderThirdCond_buf.vertL;
			in >> borderThirdCond_buf.vertR;
			in >> borderThirdCond_buf.heightL;
			in >> borderThirdCond_buf.heightR;
			in >> borderThirdCond_buf.betta;
			borderThirdCondAr[i] = borderThirdCond_buf;
		}
		return 0;

		
		return 0;
	}
	*/
};
