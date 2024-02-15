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

#pragma once

#include "visibility_solver.h"
#include "visibility_aperture_finder.h"

#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "math_plucker_2.h"
#include "math_predicates.h"
#include "plucker_polyhedron.h"
#include "plucker_polytope.h"
#include "plucker_polytope_complex.h"
#include "plucker_polytope_splitter.h"
#include "visibility_exact_query.h"
#include "silhouette.h"
#include <fstream>
#include <sstream>
#include <map>
namespace visilib
{
    template<class P>
    class  PluckerPolytope;

    /** @brief Exact visibility determination algorithm based on an aperture finder heuristic. The algorithm stops as soon as an aperture was detected.

    This class manage the CSG operations sheduling in Plucker space.
    It implements an algorithm that searches if at least on aperture exist between the two source polygons, such that at least one line that is intersecting
    the two source polygons is not blocked by any occluder. */

    template<class P, class S>
    class VisibilityApertureFinder : public VisibilitySolver<P, S>
    {
    public:
        VisibilityApertureFinder(VisibilityExactQuery_<P, S>* aSolver,
            bool normalization,
            S tolerance,
            bool detectApertureOnly);

        VisibilityResult resolve();
        VisibilityResult resolveNonRecursive();
    private:
        VisibilityResult resolveInternal(PluckerPolytope<P>* aPolytope, const std::vector<Silhouette*>& anOccluders, const std::vector<P>& aPolytopeLines, int depth
#ifdef OUTPUT_DEBUG_FILE
            , const std::string& occlusionTreeNodeSymbol
#endif
        );
        void resize(size_t myInitiaLineCount, PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P>* aPolytope);
        void extractStabbingLines(PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P>* aPolytope);
        void releasePolytope(PluckerPolytope<P>*& polytope);

        bool mNormalization;
        bool mDetectApertureOnly;
        S mTolerance;
    };

    template<class P, class S>
    VisibilityApertureFinder<P, S>::VisibilityApertureFinder(VisibilityExactQuery_<P, S>* mQuery, bool normalization, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
        mNormalization(normalization),
        mTolerance(tolerance),
        mDetectApertureOnly(detectApertureOnly)
    {
    }

    template <class P, class S>
    VisibilityResult VisibilityApertureFinder<P, S>::resolve()
    {
        return resolveInternal(reinterpret_cast<PluckerPolytope<P>*>(VisibilitySolver<P, S>::mQuery->getComplex()->getRoot()), std::vector<Silhouette*>(), std::vector<P>(), 0
#ifdef OUTPUT_DEBUG_FILE
            , "*"
#endif
        );
    }

    template<class P, class S>
    VisibilityResult VisibilityApertureFinder<P, S>::resolveInternal(PluckerPolytope<P>* aPolytope, const std::vector<Silhouette*>& anOccluders, const std::vector<P>& aPolytopeLines, int aDepth
#ifdef OUTPUT_DEBUG_FILE
        , const std::string& occlusionTreeNodeSymbol
#endif
        )
    {
        stats.maxRecursionLevel = max(aDepth, stats.maxRecursionLevel);
        stats.resolveInternal++;

        PluckerPolyhedron<P>* myPolyhedron = reinterpret_cast<PluckerPolyhedron<P>*> (VisibilitySolver<P, S>::mQuery->getComplex()->getPolyhedron());
        size_t myInitiaLineCount = myPolyhedron->getLinesCount();

#ifdef OUTPUT_DEBUG_FILE
        std::ofstream& debugOutput = VisibilitySolver<P, S>::mDebugger->getDebugOutput();
        V_LOG(debugOutput, "Resolve internal ", occlusionTreeNodeSymbol);
        aPolytope->outputProperties(debugOutput, polyhedron);
#endif
        VisibilityResult globalResult = HIDDEN;
        if (aDepth > 2000)
        {
            std::cerr << "Overflow detected: more than 2000 recursive calls...." << std::endl;
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "RESULT FAILURE: Overflow detected : more than 2000 recursive calls....", occlusionTreeNodeSymbol);
#endif
            V_ASSERT(0);
            return FAILURE;
        }
        {
            HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), STABBING_LINE_EXTRACTION);
            aPolytope->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);
        }
        if (!aPolytope->containsRealLines())
        {
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "RESULT UNKNOWN: polytope does not contains real line", occlusionTreeNodeSymbol);
#endif
            return HIDDEN;
        }

        V_ASSERT(aPolytope->isValid(myPolyhedron, mNormalization, mTolerance));

        bool hasRay = false;
        std::vector<Silhouette*> myOccluders = anOccluders;
        std::vector<P> polytopeLines = aPolytopeLines;
        if (myOccluders.empty())
        {
            if (!VisibilitySolver<P, S>::mQuery->collectAllOccluders(aPolytope, myPolyhedron, myOccluders, polytopeLines))
            {
                // Early stop - an aperture has been found: the source polygons are mutually visible
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "RAY: IS VISIBLE", occlusionTreeNodeSymbol);
#endif
                hasRay = true;
                globalResult = VISIBLE;
                if (mDetectApertureOnly)
                {
                    return VISIBLE;
                }
            }
            else
            {
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "RAY: IS OCCLUDED", occlusionTreeNodeSymbol);
#endif
            }
        }

        {
            HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
            if (VisibilitySolver<P, S>::mQuery->isOccluded(aPolytope, myPolyhedron, myOccluders, polytopeLines))
            {
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "COUNT: IS OCCLUDED", occlusionTreeNodeSymbol);
#endif
                return HIDDEN;
            }
            else
            {
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "COUNT: maybe VISIBLE", occlusionTreeNodeSymbol);
#endif
            }
        }
        Silhouette* mySilhouette = nullptr;
        size_t mySilhouetteEdgeIndex = 0;

        bool hasEdge = false;
        bool isVisible = false;

        hasEdge = VisibilitySolver<P, S>::mQuery->findNextEdge(mySilhouetteEdgeIndex, mySilhouette
            // , aPolytope
#ifdef OUTPUT_DEBUG_FILE
            , occlusionTreeNodeSymbol
#endif
        );

        if (hasEdge) // a candidate edge has been found. we will split the polytope with the hyperplane of this edge.
        {
            SilhouetteEdge& myVisibilitySilhouetteEdge = mySilhouette->getEdge(mySilhouetteEdgeIndex);
            SilhouetteMeshFace* face = myVisibilitySilhouetteEdge.mFace;
#ifndef ACTIVE_SET
            V_ASSERT(myVisibilitySilhouetteEdge.mIsActive);
#endif
            // deactivate the candidate edge for further recursion
            mySilhouette->setEdgeActive(mySilhouetteEdgeIndex, false);

            MathVector2i edge = face->getEdge(myVisibilitySilhouetteEdge.mEdgeIndex);

            // Create the Plucker representation of the edge and add hyperplane of the edge and add it to the polyhedron

            MathVector3d a = convert<MathVector3d>(face->getVertex(edge.x));
            MathVector3d b = convert<MathVector3d>(face->getVertex(edge.y));

            bool intersect = false;

            {
                HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
                intersect = MathGeometry::isEdgeInsidePolytope(a, b, aPolytope, VisibilitySolver<P, S>::mQuery->getApproximateNormal(), myPolyhedron, mTolerance);
            }
            if (intersect)
            {
                if (VisibilitySolver<P, S>::mDebugger != nullptr)
                {
                    VisibilitySolver<P, S>::mDebugger->addRemovedEdge(face->getVertex(edge.x), face->getVertex(edge.y));
                }

                size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;

                if (myPolyhedronFace == 0)
                {
                    P myHyperplane(a, b);
                    if (mNormalization)
                    {
                        myHyperplane = myHyperplane.getNormalized();
                    }
                    myPolyhedronFace = myPolyhedron->add(myHyperplane, ON_BOUNDARY, mNormalization, mTolerance);
                    myVisibilitySilhouetteEdge.mHyperPlaneIndex = myPolyhedronFace;
                }

                P myHyperplane = myPolyhedron->get(myPolyhedronFace);

                GeometryPositionType myResult = ON_NEGATIVE_SIDE;

                PluckerPolytope<P>* myPolytopeNegative = new PluckerPolytope<P>();
                PluckerPolytope<P>* myPolytopePositive = new PluckerPolytope<P>();

                {
                    HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);
                    VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SPLIT_COUNT);

#ifdef OUTPUT_DEBUG_FILE
                    V_LOG(debugOutput, "PERFORM THE SPLIT", occlusionTreeNodeSymbol);
#endif
                    myResult = PluckerPolytopeSplitter<P, S>::split(myPolyhedron, myHyperplane, aPolytope, myPolytopeNegative, myPolytopePositive, myPolyhedronFace, mNormalization, mTolerance);
#if 0
                    if (VisibilitySolver<P, S>::mQuery->getStatistic()->get(POLYTOPE_SPLIT_COUNT) % 10000 == 0)
                    {
                        VisibilitySolver<P, S>::mQuery->getStatistic()->displayCounts();
                        std::cout << std::endl;
                    }
#endif
                }

                std::vector <PluckerPolytope<P>*> myPolytopes;
#ifdef OUTPUT_DEBUG_FILE
                std::vector <std::string> postFix;
#endif
                std::vector<bool> reuseOccluders;

                if (myResult == ON_BOUNDARY)
                {
                    V_ASSERT(myPolytopeNegative && myPolytopePositive);
#ifdef OUTPUT_DEBUG_FILE
                    V_LOG(debugOutput, "SPLIT SUCCESS ->recurse on Left and right childs", occlusionTreeNodeSymbol);
#endif
                    auto position = MathPredicates::getVertexPlaneRelativePosition(myHyperplane, polytopeLines[0], mTolerance);

                    // left split
                    reuseOccluders.push_back(position != ON_POSITIVE_SIDE);
                    myPolytopes.push_back(myPolytopeNegative);
#ifdef OUTPUT_DEBUG_FILE
                    postFix.push_back("L");
#endif

                    // right split
                    reuseOccluders.push_back(position != ON_NEGATIVE_SIDE);
                    myPolytopes.push_back(myPolytopePositive);
#ifdef OUTPUT_DEBUG_FILE
                    postFix.push_back("R");
#endif
                }
                else
                {
#ifdef OUTPUT_DEBUG_FILE
                    V_LOG(debugOutput, "NO SPLIT OCCURS -> recurse again on the same polytope", occlusionTreeNodeSymbol);
#endif
                    reuseOccluders.push_back(true);
                    myPolytopes.push_back(aPolytope);
#ifdef OUTPUT_DEBUG_FILE
                    postFix.push_back("*");
#endif
                    delete myPolytopeNegative;  myPolytopeNegative = nullptr;
                    delete myPolytopePositive; myPolytopePositive = nullptr;
                }

                for (size_t i = 0; i < myPolytopes.size(); i++)
                {
                    bool hasBeenAdded = false;

                    if (myResult == ON_BOUNDARY && i == 0 || myResult == ON_NEGATIVE_SIDE)
                    {
                        mySilhouette->pushEdgeProcessed(mySilhouetteEdgeIndex);
                        hasBeenAdded = true;
                    }


#ifdef OUTPUT_DEBUG_FILE
                    std::stringstream ss;
                    ss << occlusionTreeNodeSymbol << "|" << face->getFaceIndex() << "-" << myVisibilitySilhouetteEdge.mEdgeIndex << postFix[i];
#endif
                    VisibilityResult result;

                    if (reuseOccluders[i])
                    {
                        result = resolveInternal(myPolytopes[i], myOccluders, polytopeLines, aDepth + 1
#ifdef OUTPUT_DEBUG_FILE
                            , ss.str()
#endif
                        );
                    }
                    else
                    {
                        result = resolveInternal(myPolytopes[i], std::vector<Silhouette*>(), std::vector<P>(), aDepth + 1
#ifdef OUTPUT_DEBUG_FILE
                            , ss.str()
#endif
                        );
                    }
                    if (result == FAILURE) {
                        return result;
                    }

                    if (result == VISIBLE)
                    {
                        globalResult = VISIBLE;
                        if (mDetectApertureOnly) // Early stop - an aperture has been found
                        {
#ifdef OUTPUT_DEBUG_FILE
                            V_LOG(debugOutput, "EARLY STOP - aperture found");
#endif
                            delete myPolytopeNegative;  myPolytopeNegative = nullptr;
                            delete myPolytopePositive; myPolytopePositive = nullptr;

                            return result;
                        }
                    }
                    if (hasBeenAdded)
                        mySilhouette->popEdgeProcessed(mySilhouetteEdgeIndex);

                }
                delete myPolytopeNegative;  myPolytopeNegative = nullptr;
                delete myPolytopePositive; myPolytopePositive = nullptr;

            }
            else
            {

#ifdef OUTPUT_DEBUG_FILE
                std::stringstream ss;
                ss << occlusionTreeNodeSymbol << "|" << face->getFaceIndex() << "-" << myVisibilitySilhouetteEdge.mEdgeIndex << "*";
#endif
                VisibilityResult result = resolveInternal(aPolytope, anOccluders, aPolytopeLines, aDepth + 1
#ifdef OUTPUT_DEBUG_FILE
                    , ss.str()
#endif
                );
                if (result == FAILURE)
                    return result;

                if (result == VISIBLE)
                {
                    globalResult = VISIBLE;
                    if (mDetectApertureOnly) // Early stop - an aperture has been found
                    {
                        return result;
                    }
                }
            }
            // reactivate the  edge

            mySilhouette->setEdgeActive(mySilhouetteEdgeIndex, true);
        }

        // No valid candidate edge has been found in the occluder set: all the set of lines represented by the polytope are blocked by at least one occluder: the recursion can stops.
        // When it is the case for all the sub-polytopes created during the recursive traversal, the two polygons are mutually hidden.

        if (!mDetectApertureOnly)
        {
            if (!hasEdge)
            {
                //We found a final polytope visible
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "visible polytope found: exporting stabbing lines", occlusionTreeNodeSymbol);
#endif
                extractStabbingLines(myPolyhedron, aPolytope);
                globalResult = VISIBLE;
            }
        }
#if 0
        if (polyhedron->getLinesCount() != initiaLineCount)
        {
            //     resize(initiaLineCount, polyhedron, aPolytope);
        }
#endif

        return globalResult;
    }

    template <class P, class S>
    void VisibilityApertureFinder<P, S>::releasePolytope(PluckerPolytope<P>*& polytope)
    {
        if (polytope != VisibilitySolver<P, S>::mQuery->getComplex()->getRoot())
        {
            delete polytope;
            polytope = NULL;
        }
    }

    template <class P, class S>
    VisibilityResult VisibilityApertureFinder<P, S>::resolveNonRecursive()
    {
        VisibilityResult globalResult = VisibilityResult::HIDDEN;

        PluckerPolytopeComplex<P>* complex = mQuery->getComplex();

        // get occluders
        PluckerPolyhedron<P>* polyhedron = reinterpret_cast<PluckerPolyhedron<P>*> (complex->getPolyhedron());
        size_t initiaLineCount = polyhedron->getLinesCount();

        std::vector <PluckerPolytope<P>*> polytopes;
        polytopes.push_back(reinterpret_cast<PluckerPolytope<P>*>(complex->getRoot()));

        while (!polytopes.empty())
        {
            PluckerPolytope<P>* polytope = polytopes.back();
            polytopes.pop_back();

            {
                HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), STABBING_LINE_EXTRACTION);
                polytope->computeEdgesIntersectingQuadric(polyhedron, mTolerance);
            }
            if (polytope->containsRealLines())
            {
                V_ASSERT(polytope->isValid(polyhedron, mNormalization, mTolerance));

                std::vector<Silhouette*> occluders;
                std::vector<P> polytopeLines;

                bool hasRay = false;
                if (!VisibilitySolver<P, S>::mQuery->collectAllOccluders(polytope, polyhedron, occluders, polytopeLines))
                {
                    // Early stop - an aperture has been found: the source polygons are mutually visible
                    hasRay = true;
                    globalResult = VISIBLE;
                    if (mDetectApertureOnly)
                    {
                        releasePolytope(polytope);
                        // trivially accept polytope
                        return VISIBLE;
                    }
                }
#if true
                {
                    HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
                    if (VisibilitySolver<P, S>::mQuery->isOccluded(polytope, polyhedron, occluders, polytopeLines))
                    {
                        // early out, no vertices outside of any edge of the occluder, so no subdivision required
                        releasePolytope( polytope);
                        // trivially discard polytope
                        continue;
                    }
                }
#endif

                bool visible = false;
                // iterate over all silhouettes and edges, push the unoccluded polytops onto a stack, discard the final fully occluded polytopes
                for (std::vector<Silhouette*>::const_iterator it = occluders.begin(); it != occluders.end() && polytope != nullptr; it++)
                {
                    Silhouette* mySilhouette = *it;
                    std::vector<SilhouetteEdge>& edges = mySilhouette->getEdges();
                    bool nextSilhouette = false;
                    for (size_t mySilhouetteEdgeIndex = 0; mySilhouetteEdgeIndex < edges.size() && !nextSilhouette && polytope != nullptr; mySilhouetteEdgeIndex++)
                    {
#if true
                        {
                            HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
                            if (VisibilitySolver<P, S>::mQuery->isOccluded(polytope, polyhedron, occluders, polytopeLines))
                            {
                                // early out, no vertices outside of any edge of the occluder, so no subdivision required
                                releasePolytope(polytope);
                                // trivially discard polytope
                                continue;
                            }
                        }
#endif

                        SilhouetteEdge& myVisibilitySilhouetteEdge = edges[mySilhouetteEdgeIndex];
                        SilhouetteMeshFace* face = myVisibilitySilhouetteEdge.mFace;

                        MathVector2i edge = face->getEdge(myVisibilitySilhouetteEdge.mEdgeIndex);

                        // Create the Plucker representation of the edge and add hyperplane of the edge and add it to the polyhedron

                        MathVector3d a = convert<MathVector3d>(face->getVertex(edge.x));
                        MathVector3d b = convert<MathVector3d>(face->getVertex(edge.y));

                        bool intersect = false;
#if true
                        {
                            HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
                            intersect = MathGeometry::isEdgeInsidePolytope(a, b, polytope, VisibilitySolver<P, S>::mQuery->getApproximateNormal(), polyhedron, mTolerance);
                        }
#else
                        intersect = true;
#endif
                        if (intersect)
                        {
                            if (VisibilitySolver<P, S>::mDebugger != nullptr)
                            {
                                VisibilitySolver<P, S>::mDebugger->addRemovedEdge(face->getVertex(edge.x), face->getVertex(edge.y));
                            }

                            size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;

                            if (myPolyhedronFace == 0)
                            {
                                P myHyperplane(a, b);
                                if (mNormalization)
                                {
                                    myHyperplane = myHyperplane.getNormalized();
                                }
                                myPolyhedronFace = polyhedron->add(myHyperplane, ON_BOUNDARY, mNormalization, mTolerance);
                                myVisibilitySilhouetteEdge.mHyperPlaneIndex = myPolyhedronFace;
                            }

                            P myHyperplane = polyhedron->get(myPolyhedronFace);

                            GeometryPositionType myResult = ON_NEGATIVE_SIDE;

                            PluckerPolytope<P>* myPolytopeNegative = new PluckerPolytope<P>();
                            PluckerPolytope<P>* myPolytopePositive = new PluckerPolytope<P>();

                            {
                                HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);
                                VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SPLIT_COUNT);

                                myResult = PluckerPolytopeSplitter<P, S>::split(polyhedron, myHyperplane, polytope, myPolytopeNegative, myPolytopePositive, myPolyhedronFace, mNormalization, mTolerance);
                            }

                            GeometryPositionType position = MathPredicates::getVertexPlaneRelativePosition(myHyperplane, polytopeLines[0], mTolerance);

                            switch (myResult)
                            {
                            case ON_BOUNDARY: {
                                V_ASSERT(myPolytopeNegative && myPolytopePositive);
                                if (position != ON_POSITIVE_SIDE)
                                {
                                    // we push the positive polytope into a list to process later
                                    polytopes.push_back(myPolytopePositive);
                                    stats.polytopes++;

                                    // we process the negative polytope further
                                    releasePolytope(polytope);
                                    polytope = myPolytopeNegative;
                                }
                                else
                                {
                                    // we push the positive polytope into a list to process later
                                    polytopes.push_back(myPolytopeNegative);
                                    stats.polytopes++;

                                    // we process the negative polytope further
                                    releasePolytope(polytope);
                                    polytope = myPolytopePositive;
                                }
                                break;
                            }
#if 0
                            case ON_POSITIVE_SIDE: {
                                // is on the positive side, so no point looking at the rest of the silhouetter
                                delete myPolytopeNegative; myPolytopeNegative = nullptr;
                                delete myPolytopePositive; myPolytopePositive = nullptr;
                                if (position != ON_NEGATIVE_SIDE)
                                {
                                    nextSilhouette = true;
                                    visible = true;
                                }
                                break;
                            }
                            case ON_NEGATIVE_SIDE: {
                                // is fully on the negative side of an edge, so we continue chopping bits of the existing polytope
                                delete myPolytopeNegative; myPolytopeNegative = nullptr;
                                delete myPolytopePositive; myPolytopePositive = nullptr;
                                if (position != ON_POSITIVE_SIDE)
                                {
                                    nextSilhouette = true;
                                    visible = true;
                                }
                                break;
                            }
#else
                            default: {
                                delete myPolytopeNegative; myPolytopeNegative = nullptr;
                                delete myPolytopePositive; myPolytopePositive = nullptr;
                                break;
                            }
#endif
                            }
                        }
                    }
                }
                if (!mDetectApertureOnly && visible)
                {
                    // this produces the green lines:
                    polytope->computeEdgesIntersectingQuadric(polyhedron, mTolerance);
                    extractStabbingLines(polyhedron, polytope);
                }
            }

            releasePolytope(polytope);
            polytope = nullptr;
        }
        return globalResult;
    }

    template<class P, class S>
    void VisibilityApertureFinder<P, S>::resize(size_t anInitiaLineCount, PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P>* aPolytope)
    {
        HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);

        myPolyhedron->resize(anInitiaLineCount);
        for (auto v : aPolytope->getVertices())
        {
            std::vector<size_t>& facets = myPolyhedron->getFacetsDescription(v);
            for (size_t i = 0; i < facets.size(); i++)
            {
                if (facets[i] >= anInitiaLineCount)
                {
                    facets.resize(i);
                }
            }
        }
    }

    template<class P, class S>
    void VisibilityApertureFinder<P, S>::extractStabbingLines(PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P>* aPolytope)
    {
        HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), STABBING_LINE_EXTRACTION);

        MathPlane3d aPlane0 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(0)->getPlane();
        MathPlane3d aPlane1 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(1)->getPlane();

        if (aPolytope->getExtremalStabbingLinesCount() == 0)
        {
            aPolytope->computeExtremalStabbingLines(myPolyhedron, mTolerance);
        }

        if (VisibilitySolver<P, S>::mDebugger != nullptr)
        {
            std::vector<std::pair<MathVector3d, MathVector3d>> lines;
            aPolytope->getExtremalStabbingLinesBackTo3D(lines, aPlane0, aPlane1);
            //for (auto line : lines)
            for (size_t lineIndex = 0; lineIndex < lines.size(); lineIndex++)
            {
                const std::pair<MathVector3d, MathVector3d>& line = lines[lineIndex];
                VisibilitySolver<P, S>::mDebugger->addExtremalStabbingLine(convert<MathVector3f>(line.first), convert<MathVector3f>(line.second));
            }
        }
    }
}
