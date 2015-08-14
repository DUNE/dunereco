#ifndef IniSegAlg_h
#define IniSegAlg_h

#include "RecoAlg/ProjectionMatchingAlg.h"

namespace fhicl {
	class ParameterSet;
}

namespace dunefd {
	class Hit2D;
	class IniSegAlg;
}

class dunefd::Hit2D
{
	public:
	Hit2D(TVector2 point2d, size_t key);

	private:
	TVector2 fPoint;
	size_t fKey;
};

class dunefd::IniSegAlg
{
	public:
	IniSegAlg();

	private:

};

#endif //IniSegAlg_h
