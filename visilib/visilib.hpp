/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2021 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/


#include "visilib_core.h"
#include "visilib.h"
#include "visibility_exact_query.h"
#include "geometry_convex_polygon.h"
#include "geometry_occluder_set.h"

#ifdef ENABLE_GMP
#include <gmp.h>
#endif

using namespace visilib;

inline VisibilityResult visilib::areVisible(GeometryOccluderSet* scene, const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1,
    const VisibilityExactQueryConfiguration& configuration, HelperVisualDebugger* debugger)
{
    if (vertices0 == 0 || vertices1 == 0)
    {
        std::cerr << "Error: invalid number of vertices" << std::endl;
        return FAILURE;
    }

    if (scene == nullptr || dynamic_cast<GeometryOccluderSet*>(scene) == nullptr)
    {
        std::cerr << "Error: invalid scene" << std::endl;
        return FAILURE;
    }
    if (vertices0 == nullptr || vertices1 == nullptr)
    {
        std::cerr << "Error: invalid vertex array" << std::endl;
        return FAILURE;
    }

    VisibilityExactQuery* query = nullptr;

    switch (configuration.precision)
    {
#ifdef EXACT_ARITHMETIC
    case VisibilityExactQueryConfiguration::EXACT:
        query = new VisibilityExactQuery_<MathPlucker6<exact>, exact>(scene, configuration, MathArithmetic<exact>::Tolerance());
        break;
#endif
#ifdef ENABLE_GMP
    case VisibilityExactQueryConfiguration::GMP_FLOAT:
        query = new VisibilityExactQuery_<MathPlucker6<GmpFloat>, GmpFloat>(scene, configuration, MathArithmetic<GmpFloat>::Tolerance());
        break;
    case VisibilityExactQueryConfiguration::GMP_RATIONAL:
        query = new VisibilityExactQuery_<MathPlucker6<GmpRational>, GmpRational>(scene, configuration, MathArithmetic<GmpRational>::Tolerance());
        break;
#endif
#ifdef ENABLE_REALEXPR
    case VisibilityExactQueryConfiguration::REAL_EXPR:
        query = new VisibilityExactQuery_<MathPlucker6<RealExpr>, RealExpr>(scene, configuration, MathArithmetic<RealExpr>::Tolerance());
        break;
#endif
#ifdef ENABLE_MPFR
    case VisibilityExactQueryConfiguration::MPFR:
        query = new VisibilityExactQuery_<MathPlucker6<Mpfr>, Mpfr>(scene, configuration, MathArithmetic<Mpfr>::Tolerance());
        break;
#endif
    case VisibilityExactQueryConfiguration::DOUBLE:
        query = new VisibilityExactQuery_<MathPlucker6<double>, double>(scene, configuration, MathArithmetic<double>::Tolerance());
        break;

    default:

        query = new VisibilityExactQuery_<MathPlucker6<float>, float>(scene, configuration, MathArithmetic<float>::Tolerance());
        break;
    }

    query->attachVisualisationDebugger(debugger);

    VisibilityResult result = query->arePolygonsVisible(vertices0, numVertices0, vertices1, numVertices1);

    if (debugger)
        query->displayStatistic();
    delete query;

    return result;
}

