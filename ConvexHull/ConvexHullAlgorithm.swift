//
//  ConvexHullAlgorithm.swift
//  ConvexHull
//
//  Created by Carl Wieland on 1/18/18.
/******************************************************************************
 *
 * The MIT License (MIT)
 *
 * MIConvexHull, Copyright (c) 2015 David Sehnal, Matthew Campbell
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *****************************************************************************/

import Foundation

internal class ConvexHullAlgorithm {

    /// Are we on a paraboloid?
    private let isLifted: Bool

    /// Explained in ConvexHullComputationConfig.
    private let planeDistanceTolerance: Double

    /*
     * Representation of the input vertices.
     *
     * - In the algorithm, a vertex is represented by its index in the Vertices array.
     *   This makes the algorithm a lot faster (up to 30%) than using object reference everywhere.
     * - Positions are stored as a single array of values. Coordinates for vertex with index i
     *   are stored at indices <i * Dimension, (i + 1) * Dimension)
     * - VertexMarks are used by the algorithm to help identify a set of vertices that is "above" (or "beyond")
     *   a specific face.
     */
    /// The vertices
    private let vertices: [Vector3]

    /// The positions
    private var positions: SimpleList<Double>

    /// The vertex marks
    private var vertexVisited: [Bool]

    private let numberOfVertices: Int

    /*
     * The triangulation faces are represented in a single pool for objects that are being reused.
     * This allows for represent the faces as integers and significantly speeds up many computations.
     * - AffectedFaceFlags are used to mark affected faces/
     */
    /// The face pool
    internal let facePool = SimpleList<ConvexFaceInternal>()

    /// The affected face flags
    internal var affectedFaceFlags: [Bool]

    /// Used to track the size of the current hull in the Update/RollbackCenter functions.
    private var convexHullSize = 0

    /// A list of faces that that are not a part of the final convex hull and still need to be processed.
    private let unprocessedFaces = FaceList()

    /// A list of faces that form the convex hull.
    private let convexFaces = IndexBuffer()

    /// The vertex that is currently being processed.
    private var currentVertex = 0

    /// A helper variable to determine the furthest vertex for a particular convex face.
    private var  maxDistance = 0.0

    /// A helper variable to help determine the index of the vertex that is furthest from the face that is currently being
    /// processed.
    private var furthestVertex = 0


    /// The centroid of the currently computed hull.
    private var center = Vector3.zero

    /*
     * Helper arrays to store faces for adjacency update.
     * This is just to prevent unnecessary allocations.
     */

    /// The update buffer
    private var updateBuffer: [Int]

    /// The update indices
    private var updateIndices: [Int]


    /// Used to determine which faces need to be updated at each step of the algorithm.
    private let traverseStack = IndexBuffer()


    /// Used for VerticesBeyond for faces that are on the convex hull.
    private let emptyBuffer = IndexBuffer()


    /// Used to determine which vertices are "above" (or "beyond") a face
    private var beyondBuffer = [Int]()


    /// Stores faces that are visible from the current vertex.
    private let affectedFaceBuffer = IndexBuffer()


    /// Stores faces that form a "cone" created by adding new vertex.
    private let coneFaceBuffer = SimpleList<DeferredFace>()

    /// Stores a list of "singular" (or "generate", "planar", etc.) vertices that cannot be part of the hull.
    private var singularVertices = Set<Int>()


    /// The connector table helps to determine the adjacency of convex faces.
    /// Hashing is used instead of pairwise comparison. This significantly speeds up the computations,
    /// especially for higher dimensions.
    private var connectorTable = [UInt64: [FaceConnector]]()


    /// Manages the memory allocations and storage of unused objects.
    private let objectManager: ObjectManager

    /// Helper class for handling math related stuff.
    private let mathHelper: MathHelper
    private var boundingBoxPoints: [[Int]]
    private var indexOfDimensionWithLeastExtremes = 0
    private var minima = [0.0,0.0,0.0]
    private var maxima = [0.0,0.0,0.0]



    public static func getConvexHull(with data:[Vector3], planeDistanceTolerance tolerance: Double) -> ConvexHull {

        let ch = ConvexHullAlgorithm(vertices: data, lift: false, planeDistanceTolerance: tolerance)
        ch.generateConvexHull()
        return ConvexHull(points: ch.hullVertices(data: data), faces: ch.getConvexFaces())
    }

    private init( vertices: [Vector3],  lift: Bool, planeDistanceTolerance: Double) {
        isLifted = lift
        self.vertices = vertices
        numberOfVertices = vertices.count


        boundingBoxPoints = [[Int]](repeating: [], count: 3)
        vertexVisited = [Bool](repeating: false, count: numberOfVertices)
        positions = SimpleList<Double>()
        affectedFaceFlags = [Bool](repeating: false, count: (3 + 1) * 10)
        updateBuffer = [Int](repeating: 0, count: 3)
        updateIndices = [Int](repeating: 0, count: 3)

        self.planeDistanceTolerance = planeDistanceTolerance

        mathHelper = MathHelper(positions: positions);
        objectManager = ObjectManager(facePool: facePool)

        repeatElement(0, count: numberOfVertices * 3).forEach({ positions.append($0)})

    }

    deinit {
        
    }

    /// Gets/calculates the convex hull. This is
    private func generateConvexHull() {
        // accessing a 1D array is quicker than a jagged array, so the first step is to make this array
        serializeVerticesToPositions();
        // next the bounding box extremes are found. This is used to shift, scale and find the starting simplex.
        findBoundingBoxPoints();
        // the positions are shifted to avoid divide by zero problems
        // and if Delaunay or Voronoi, then the parabola terms are scaled back to match the size of the other coords
        shiftAndScalePositions();
        // Find the (dimension+1) initial points and create the simplexes.
        createInitialSimplex();

        // Now, the main loop. These initial faces of a simplex are replaced and expanded
        // outwards to make the convex hull and faces.
        while let currentFace = unprocessedFaces.first {

            currentVertex = currentFace.furthestVertex;

            updateCenter();

            // The affected faces get tagged
            tagAffectedFaces(from: currentFace);

            // Create the cone from the currentVertex and the affected faces horizon.
            if !singularVertices.contains(currentVertex) && createCone() {
                commitCone()
            } else {
                handleSingular();
            }

            // Need to reset the tags
            let count = affectedFaceBuffer.count;
            for i in 0..<count {
                affectedFaceFlags[affectedFaceBuffer[i]] = false
            }
        }
    }

    /// Serializes the vertices into the 1D array, Positions. The 1D array has much quicker access
    private func serializeVerticesToPositions() {
        var index = 0;
        if isLifted { // "Lifted" means that the last dimension is the sum of the squares of the others.
            for v in vertices {
                var parabolaTerm = 0.0 // the lifted term is a sum of squares.
                let origNumDim = 3 - 1;
                for i in 0..<origNumDim {
                    let coordinate = v.position[i]
                    positions[index] = coordinate;
                    index += 1
                    parabolaTerm += coordinate * coordinate;
                }
                positions[index] = parabolaTerm;
                index += 1
            }
        } else {
            for v in vertices {
                for i in 0..<3 {
                    positions[index] = v.position[i]
                    index += 1
                }
            }
        }
    }

    /// Finds the bounding box points.
    private func findBoundingBoxPoints() {
        indexOfDimensionWithLeastExtremes = -1;
        var minNumExtremes = Int.max
        for i in 0..<3 {
            var minIndices = [Int]()
            var maxIndices = [Int]()
            var min = Double.greatestFiniteMagnitude
            var max = Double.greatestFiniteMagnitude * -1

            for j in 0..<numberOfVertices {
                let v = getCoordinate(vIndex: j, dimension: i);
                var difference = min - v;
                if difference >= planeDistanceTolerance {
                    // you found a better solution than before, clear out the list and store new value
                    min = v;
                    minIndices.removeAll(keepingCapacity: true)
                    minIndices.append(j);
                } else if (difference > 0) {
                    // you found a solution slightly better than before, clear out those that are no longer on the list and store new value
                    min = v;
                    minIndices.removeWhere { min - getCoordinate(vIndex: $0, dimension: i) > planeDistanceTolerance }
                    minIndices.append(j);
                } else if (difference > -planeDistanceTolerance) {
                    //same or almost as good as current limit, so store it
                    minIndices.append(j);
                }

                difference = v - max;

                if difference >= planeDistanceTolerance {
                    // you found a better solution than before, clear out the list and store new value
                    max = v;
                    maxIndices.removeAll(keepingCapacity: true)
                    maxIndices.append(j);
                } else if difference > 0 {
                    // you found a solution slightly better than before, clear out those that are no longer on the list and store new value
                    max = v;
                    maxIndices.removeWhere { min - getCoordinate(vIndex: $0, dimension: i) > planeDistanceTolerance }
                    maxIndices.append(j);
                } else if difference > -planeDistanceTolerance {
                    //same or almost as good as current limit, so store it
                    maxIndices.append(j);
                }
            }
            minima[i] = min;
            maxima[i] = max;
            minIndices.append(contentsOf: maxIndices)
            if minIndices.count < minNumExtremes {
                minNumExtremes = minIndices.count;
                indexOfDimensionWithLeastExtremes = i;
            }
            boundingBoxPoints[i] = minIndices;
        }
    }



    /// Shifts and scales the Positions to avoid future errors. This does not alter the original data.
    private func shiftAndScalePositions() {
        let positionsLength = positions.count
        if isLifted {
            let origNumDim = 3 - 1;
            let minSum = Double(minima.reduce(0, {$0 + abs($1)}))
            let maxSum = Double(maxima.reduce(0, { $0 + abs($1)}))
            let parabolaScale = 2.0 / (minSum + maxSum - abs(maxima[origNumDim]) - abs(minima[origNumDim]))
            // the parabola scale is 1 / average of the sum of the other dimensions.
            // multiplying this by the parabola will scale it back to be on near similar size to the
            // other dimensions. Without this, the term is much larger than the others, which causes
            // problems for roundoff error and finding the normal of faces.
            minima[origNumDim] *= parabolaScale; // change the extreme values as well
            maxima[origNumDim] *= parabolaScale;
            // it is done here because
            for i in stride(from: origNumDim, to:positionsLength, by: 3) {
                positions[i] *= parabolaScale;
            }

        }
        var shiftAmount = Vector3.zero
        for i in 0..<3 {
            // now the entire model is shifted to all positive numbers...plus some more.
            // why?
            // 1) to avoid dealing with a point at the origin {0,0,...,0} which causes problems
            //    for future normal finding
            // 2) note that weird shift that is used (max - min - min). This is to avoid scaling
            //    issues. this shift means that the minima in a dimension will always be a positive
            //    number (no points at zero), and the minima [in a given dimension] will always be
            //    half of the maxima. 'Half' is much preferred to 'thousands of times'
            //    Think of the first term as the range (max - min), then the second term avoids cases
            //    where there are both positive and negative numbers.
            if (maxima[i] == minima[i]){
                shiftAmount[i] = 0.0;
            } else {
                shiftAmount[i] = (maxima[i] - minima[i]) - minima[i];
            }
        }
        for i in 0..<positionsLength {
            positions[i] += shiftAmount[i % 3];
        }

    }

    /// Find the (dimension+1) initial points and create the simplexes.
    /// Creates the initial simplex of n+1 vertices by using points from the bounding box.
    /// Special care is taken to ensure that the vertices chosen do not result in a degenerate shape
    /// where vertices are collinear (co-planar, etc). This would technically be resolved when additional
    /// vertices are checked in the main loop, but: 1) a degenerate simplex would not eliminate any other
    /// vertices (thus no savings there), 2) the creation of the face normal is prone to error.
    private func createInitialSimplex() {
        var initialPoints = findInitialPoints();
        //create the first faces from (dimension + 1) vertices.
        var faces = [Int](repeating:0, count: 3 + 1)

        for i in 0..<(3 + 1) {
            var vertices = [0,0,0]
            var k = 0
            for j in 0...3 {
                if i != j {
                    vertices[k] = initialPoints[j];
                    k += 1
                }
            }
            let newFace = facePool[objectManager.getFace()]

            newFace.vertices = vertices.sorted()

            _ = mathHelper.calculateFacePlane(face: newFace, center: center);
            faces[i] = newFace.index;
        }
        // update the adjacency (check all pairs of faces)
        for i in 0..<3 {
            for j in (i + 1)..<(3 + 1) {
                updateAdjacency(l: facePool[faces[i]], r: facePool[faces[j]]);
            }
        }


        // Init the vertex beyond buffers.

        for faceIndex in faces  {
            let face = facePool[faceIndex];
            findBeyondVertices(of: face);
            if face.verticesBeyond.count == 0 {
                convexFaces.append(face.index); // The face is on the hull
            } else {
                unprocessedFaces.append(face);
            }
        }


        // Set all vertices to false (unvisited).
        for vertex in initialPoints {
            vertexVisited[vertex] = false;
        }
    }



    /// Finds (dimension + 1) initial points.
    private func findInitialPoints() -> [Int] {
        let bigNumber = maxima.reduce(0,+) * Double(3 * numberOfVertices)
        // the first two points are taken from the dimension that had the fewest extremes
        // well, in most cases there will only be 2 in all dimensions: one min and one max
        // but a lot of engineering part shapes are nice and square and can have hundreds of
        // parallel vertices at the extremes
        let vertex1 = boundingBoxPoints[indexOfDimensionWithLeastExtremes].removeFirst() // these are min and max vertices along
        let vertex2 = boundingBoxPoints[indexOfDimensionWithLeastExtremes].removeLast() // the dimension that had the fewest points
        var initialPoints = [vertex1, vertex2 ]
        vertexVisited[vertex1] = true
        vertexVisited[vertex2] = true;
        currentVertex = vertex1
        updateCenter()
        currentVertex = vertex2
        updateCenter()
        var edgeVectors = [Vector3](repeating: .zero, count:3)
        edgeVectors[0] = mathHelper.vectorBetweenVertices(toIndex: vertex2, fromIndex: vertex1);
        // now the remaining vertices are just combined in one big list
        var extremes = boundingBoxPoints.flatMap({ $0})
        // otherwise find the remaining points by maximizing the initial simplex volume
        var index = 1;
        while index < 3 && !extremes.isEmpty {
            var bestVertex = -1;
            var bestEdgeVector = Vector3.zero
            var maxVolume = Constants.defaultPlaneDistanceTolerance;
            for i in stride(from: extremes.count - 1, through: 0, by: -1) {
                // count backwards in order to remove potential duplicates
                let vIndex = extremes[i];
                if initialPoints.contains(vIndex){
                    extremes.remove(at: i)
                } else {
                    edgeVectors[index] = mathHelper.vectorBetweenVertices(toIndex: vIndex, fromIndex: vertex1);
                    let volume = mathHelper.getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber);
                    if maxVolume < volume {
                        maxVolume = volume;
                        bestVertex = vIndex;
                        bestEdgeVector = edgeVectors[index];
                    }
                }
            }
            if let index = extremes.index(of: bestVertex) {
                extremes.remove(at: index)
            }
            if bestVertex == -1 {
                break;
            }
            initialPoints.append(bestVertex);
            edgeVectors[index] = bestEdgeVector;
            index += 1
            currentVertex = bestVertex
            updateCenter();
        }
        // hmm, there are not enough points on the bounding box to make a simplex. It is rare but entirely possibly.
        // As an extreme, the bounding box can be made in n dimensions from only 2 unique points. When we can't find
        // enough unique points, we start again with ALL the vertices. The following is a near replica of the code
        // above, but instead of extremes, we consider "allVertices".
        if initialPoints.count <= 3 && !isLifted {
            var allVertices = [Int](0..<numberOfVertices)
            while index < 3 && !allVertices.isEmpty {
                var bestVertex = -1;
                var bestEdgeVector = Vector3.zero
                var maxVolume = 0.0;
                for i in stride(from: allVertices.count - 1, through: 0, by: -1) {
                    // count backwards in order to remove potential duplicates
                    let vIndex = allVertices[i];
                    if initialPoints.contains(vIndex) {
                        allVertices.remove(at: i)
                    } else {
                        edgeVectors[index] = mathHelper.vectorBetweenVertices(toIndex: vIndex, fromIndex: vertex1);
                        let volume = mathHelper.getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber);
                        if maxVolume < volume {
                            maxVolume = volume
                            bestVertex = vIndex
                            bestEdgeVector = edgeVectors[index]
                        }
                    }
                }
                if let index = allVertices.index(of: bestVertex) {
                    allVertices.remove(at: index)
                }
                if (bestVertex == -1) {
                    break;
                }
                initialPoints.append(bestVertex)
                edgeVectors[index] = bestEdgeVector
                index += 1
                currentVertex = bestVertex
                updateCenter();
            }
        }
        if initialPoints.count <= 3 && isLifted {
            var allVertices = [Int](0..<numberOfVertices)
            while index < 3 && !allVertices.isEmpty {
                var bestVertex = -1;
                var bestEdgeVector = Vector3.zero
                var maxVolume = 0.0;
                for i in stride(from: allVertices.count - 1, through: 0, by: -1) {
                    // count backwards in order to remove potential duplicates
                    let vIndex = allVertices[i];
                    if initialPoints.contains(vIndex) {
                        allVertices.remove(at: i)
                    } else {
                        mathHelper.randomOffsetToLift(index: vIndex, maxHeight: maxima.last! - minima.last!)
                        edgeVectors[index] = mathHelper.vectorBetweenVertices(toIndex: vIndex, fromIndex: vertex1);
                        let volume = mathHelper.getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber);
                        if (maxVolume < volume)
                        {
                            maxVolume = volume;
                            bestVertex = vIndex;
                            bestEdgeVector = edgeVectors[index];
                        }
                    }
                }
                if let index = allVertices.index(of: bestVertex) {
                    allVertices.remove(at: index)
                }
                if (bestVertex == -1) {
                    break;
                }
                initialPoints.append(bestVertex)
                edgeVectors[index] = bestEdgeVector
                index += 1
                currentVertex = bestVertex
                updateCenter();
            }
        }
        if initialPoints.count <= 3 && isLifted {
            fatalError("The input data is degenerate. It appears to exist in \(3) dimensions, but it is a \(3 - 1) dimensional set (i.e. the points are collinear, coplanar, or co-hyperplanar.)")
        }
        return initialPoints;
    }


    /// Check if 2 faces are adjacent and if so, update their AdjacentFaces array.
    private func updateAdjacency(l: ConvexFaceInternal, r: ConvexFaceInternal) {
        let lv = l.vertices;
        let rv = r.vertices;

        // reset marks on the 1st face
        lv.forEach { vertexVisited[$0] = false}

        // mark all vertices on the 2nd face
        rv.forEach({ vertexVisited[$0] = true })

        var i = 0

        // find the 1st false index
        for k in 0..<lv.count {
            if !vertexVisited[lv[k]] {
                i = k
                break;
            }
        }
        // no vertex was marked
        if i == 3 {
            return;
        }


        // check if only 1 vertex wasn't marked
        for j in (i + 1)..<lv.count {
            if !vertexVisited[lv[j]] {
                return;
            }
        }

        // if we are here, the two faces share an edge
        l.adjacentFaces[i] = r.index;

        // update the adj. face on the other face - find the vertex that remains marked
        lv.forEach { vertexVisited[$0] = false}
        for k in 0..<rv.count {
            if vertexVisited[rv[k]] {
                i = k
                break
            }
        }
        r.adjacentFaces[i] = l.index;
    }


    /// Used in the "initialization" code.
    private func findBeyondVertices(of face: ConvexFaceInternal) {
        maxDistance = -Double.greatestFiniteMagnitude
        furthestVertex = 0
        for i in 0..<numberOfVertices where !vertexVisited[i] {
            isBeyond(face: face, beyondVertices: &face.verticesBeyond, v: i)
        }

        face.furthestVertex = furthestVertex;
    }

    /// Get a vertex coordinate. In order to reduce speed, all vertex coordinates
    /// have been placed in a single array.
    private func getCoordinate(vIndex: Int, dimension: Int) -> Double {
        return positions[(vIndex * 3) + dimension]
    }


    /// Check whether the vertex v is beyond the given face. If so, add it to beyondVertices.
    private func isBeyond(face: ConvexFaceInternal, beyondVertices: inout [Int], v: Int) {
        let distance = mathHelper.getVertexDistance(v: v, f: face);
        if distance >= planeDistanceTolerance {
            if distance > maxDistance {
                // If it's within the tolerance distance, use the lex. larger point
                if distance - maxDistance < planeDistanceTolerance { // todo: why is this LexCompare necessary. Would seem to favor x over y over z (etc.)?
                    if lexCompare(u: v, v: furthestVertex) > 0 {
                        maxDistance = distance;
                        furthestVertex = v;
                    }
                } else {
                    maxDistance = distance;
                    furthestVertex = v;
                }
            }
            beyondVertices.append(v);
        }
    }

    /// Compares the values of two vertices. The return value (-1, 0 or +1) are found
    /// by first checking the first coordinate and then progressing through the rest.
    /// In this way {2, 8} will be a "-1" (less than) {3, 1}.
    private func lexCompare(u: Int, v: Int) -> Int {
        let uOffset = u * 3, vOffset = v * 3
        for i in 0..<3 {
            let x = positions[uOffset + i], y = positions[vOffset + i];
            if x < y {
                return -1
            } else if x > y {
                return 1
            }
        }
        return 0;
    }

    /// Recalculates the centroid of the current hull.
    private func updateCenter() {
        for i in 0..<3 {
            center[i] *= Double(convexHullSize)
        }
        convexHullSize += 1;
        let f = 1.0 / Double(convexHullSize)
        let co = currentVertex * 3;
        for i in 0..<3 {
            center[i] = f * (center[i] + positions[co + i])
        }
    }


    /// Tags all faces seen from the current vertex with 1.
    private func tagAffectedFaces(from currentFace: ConvexFaceInternal) {
        affectedFaceBuffer.clear()
        affectedFaceBuffer.append(currentFace.index);
        traverseAffectedFaces(from: currentFace.index);
    }


    /// Recursively traverse all the relevant faces.
    private func traverseAffectedFaces(from currentFace: Int) {
        traverseStack.clear();
        traverseStack.append(currentFace)
        affectedFaceFlags[currentFace] = true
        while let toVisit = traverseStack.pop() {
            let top = facePool[toVisit]
            for i in 0..<3 {
                let adjFace = top.adjacentFaces[i];

                if !affectedFaceFlags[adjFace] && mathHelper.getVertexDistance(v: currentVertex, f: facePool[adjFace]) >= planeDistanceTolerance {
                    affectedFaceBuffer.append(adjFace);
                    affectedFaceFlags[adjFace] = true;
                    traverseStack.append(adjFace);
                }
            }
        }
    }

    /// Connect faces using a connector.
    private func connectFace(with connector: FaceConnector) {
        let index = connector.hashCode % UInt64(Constants.connectorTableSize)
        var list = connectorTable[index] ?? []
        defer {
            connectorTable[index] = list
        }

        for (index, current) in list.enumerated() {
            if FaceConnector.areConnectable(a: connector, b: current) {
                list.remove(at: index)
                FaceConnector.connect(a: current, b: connector)
                current.face = nil
                connector.face = nil
                objectManager.depositConnector(current);
                objectManager.depositConnector(connector);
                return;
            }
        }

        list.append(connector);
    }

    /// Removes the faces "covered" by the current vertex and adds the newly created ones.
    private func createCone() -> Bool {
        let currentVertexIndex = currentVertex;
        coneFaceBuffer.clear();

        for fIndex in 0..<affectedFaceBuffer.count {
            let oldFaceIndex = affectedFaceBuffer[fIndex];
            let oldFace = facePool[oldFaceIndex];

            // Find the faces that need to be updated
            var updateCount = 0;
            for i in 0..<3 {
                let af = oldFace.adjacentFaces[i]
                if !affectedFaceFlags[af] {
                    updateBuffer[updateCount] = af
                    updateIndices[updateCount] = i
                    updateCount += 1
                }
            }

            for i in 0..<updateCount {
                let adjacentFace = facePool[updateBuffer[i]];

                var oldFaceAdjacentIndex = 0;
                let adjFaceAdjacency = adjacentFace.adjacentFaces
                for j in 0..<adjFaceAdjacency.count {
                    if (oldFaceIndex == adjFaceAdjacency[j]) {
                        oldFaceAdjacentIndex = j;
                        break;
                    }
                }

                let forbidden = updateIndices[i]; // Index of the face that corresponds to this adjacent face

                var oldVertexIndex = 0

                let newFaceIndex = objectManager.getFace();
                let newFace = facePool[newFaceIndex];

                for j in 0..<3 {
                    newFace.vertices[j] = oldFace.vertices[j];
                }
                oldVertexIndex = newFace.vertices[forbidden];

                var orderedPivotIndex = 0

                // correct the ordering
                if currentVertexIndex < oldVertexIndex {
                    orderedPivotIndex = 0;
                    for j in stride(from: forbidden - 1,through: 0, by: -1) {
                        if newFace.vertices[j] > currentVertexIndex {
                            newFace.vertices[j + 1] = newFace.vertices[j];
                        } else {
                            orderedPivotIndex = j + 1;
                            break;
                        }
                    }
                } else {
                    orderedPivotIndex = 3 - 1;
                    for j in (forbidden + 1)..<3 {
                        if newFace.vertices[j] < currentVertexIndex {
                            newFace.vertices[j - 1] = newFace.vertices[j];
                        } else {
                            orderedPivotIndex = j - 1;
                            break;
                        }
                    }
                }

                newFace.vertices[orderedPivotIndex] = currentVertex;

                if !mathHelper.calculateFacePlane(face: newFace, center: center) {
                    return false;
                }
                let deferredFace = DeferredFace(face: newFace, pivot: adjacentFace, oldFace: oldFace, faceIndex: orderedPivotIndex, pivotIndex: oldFaceAdjacentIndex)

                coneFaceBuffer.append(deferredFace)

            }
        }

        return true;
    }

    /// Commits a cone and adds a vertex to the convex hull.
    private func commitCone() {
        // Fill the adjacency.
        for  i in 0..<coneFaceBuffer.count {
            let face = coneFaceBuffer[i];

            let newFace = face.face
            let adjacentFace = face.pivot
            let oldFace = face.oldFace

            let orderedPivotIndex = face.faceIndex;

            newFace.adjacentFaces[orderedPivotIndex] = adjacentFace.index;
            adjacentFace.adjacentFaces[face.pivotIndex] = newFace.index;

            // let there be a connection.
            for j in 0..<3 {
                if (j == orderedPivotIndex) {
                    continue;

                }
                let connector = objectManager.getConnector();
                connector.update(face: newFace, edgeIndex: j)
                connectFace(with: connector);
            }

            // the id adjacent face on the hull? If so, we can use simple method to find beyond vertices.
            if (adjacentFace.verticesBeyond.count == 0){
                findBeyondVertices(face: newFace, beyond: &oldFace.verticesBeyond);
            } else if (adjacentFace.verticesBeyond.count < oldFace.verticesBeyond.count) { // it is slightly more effective if the face with the lower number of beyond vertices comes first.
                findBeyondVertices(face: newFace, beyond: &adjacentFace.verticesBeyond, beyond1: &oldFace.verticesBeyond);
            } else {
                findBeyondVertices(face: newFace, beyond: &oldFace.verticesBeyond, beyond1: &adjacentFace.verticesBeyond);
            }

            // This face will definitely lie on the hull
            if newFace.verticesBeyond.isEmpty {
                convexFaces.append(newFace.index);
                unprocessedFaces.remove(newFace);
                newFace.verticesBeyond.removeAll(keepingCapacity: true)
            } else {
                unprocessedFaces.append(newFace);
            }

        }

        // Recycle the affected faces.
        for fIndex in 0..<affectedFaceBuffer.count {
            let faceIndex = affectedFaceBuffer[fIndex];
            unprocessedFaces.remove(facePool[faceIndex]);
            objectManager.depositFace(at: faceIndex);
        }
    }

    /// Handles singular vertex.
    private func handleSingular() {
        rollbackCenter();
        singularVertices.insert(currentVertex);

        // This means that all the affected faces must be on the hull and that all their "vertices beyond" are singular.
        for fIndex in 0..<affectedFaceBuffer.count {
            let face = facePool[affectedFaceBuffer[fIndex]]
            let vb = face.verticesBeyond;
            for i in 0..<vb.count {
                singularVertices.insert(vb[i])
            }

            convexFaces.append(face.index);
            unprocessedFaces.remove(face);
            face.verticesBeyond.removeAll(keepingCapacity: true)
        }
    }

    /// Removes the last vertex from the center.
    private func rollbackCenter() {
        for i in 0..<3 {
            center[i] *= Double(convexHullSize)
        }
        convexHullSize -= 1;
        let f = convexHullSize > 0 ? 1.0 / Double(convexHullSize) : 0.0
        let co = currentVertex * 3;
        for i in 0..<3 {
            center[i] = f * (center[i] - positions[co + i])
        }
    }

    /// Used by update faces.
    private func findBeyondVertices(face: ConvexFaceInternal, beyond: inout [Int] , beyond1: inout [Int] ) {
        maxDistance = Double.greatestFiniteMagnitude * -1
        furthestVertex = 0;
        var v = 0

        for i in 0..<beyond1.count {
            vertexVisited[beyond1[i]] = true;
        }

        vertexVisited[currentVertex] = false;
        for i in 0..<beyond.count {
            v = beyond[i];
            if v == currentVertex{
                continue
            }
            vertexVisited[v] = false;
            isBeyond(face: face, beyondVertices: &beyondBuffer, v: v)
        }

        for i in 0..<beyond1.count {
            v = beyond1[i];
            if vertexVisited[v] {
                isBeyond(face: face, beyondVertices: &beyondBuffer, v: v);
            }
        }

        face.furthestVertex = furthestVertex;

        // Pull the old switch a roo (switch the face beyond buffers)
        var temp = face.verticesBeyond;
        face.verticesBeyond = beyondBuffer
        if (temp.count > 0){
            temp.removeAll(keepingCapacity: true)
        }
        beyondBuffer = temp;
    }

    /// Finds the beyond vertices.
    private func findBeyondVertices(face: ConvexFaceInternal , beyond: inout [Int]){

        maxDistance = -Double.greatestFiniteMagnitude
        furthestVertex = 0;
        var v = 0

        for i in 0..<beyond.count {
            v = beyond[i];
            if v == currentVertex {
                continue
            }
            isBeyond(face: face, beyondVertices: &beyondBuffer, v: v);
        }

        face.furthestVertex = furthestVertex;

        // Pull the old switch a roo (switch the face beyond buffers)
        var temp = face.verticesBeyond
        face.verticesBeyond = beyondBuffer
        if temp.count > 0 {
            temp.removeAll(keepingCapacity: true)
        }
        beyondBuffer = temp
    }


    /// Gets the hull vertices.
    private func hullVertices(data: [Vector3]) -> [Vector3] {
        let cellCount = convexFaces.count

        for i in 0..<numberOfVertices {
            vertexVisited[i] = false
        }

        for i in 0..<cellCount {
            let vs = facePool[convexFaces[i]].vertices
            for j in 0..<vs.count {
                let v = vs[j]
                if !vertexVisited[v] {
                    vertexVisited[v] = true
                }
            }
        }

        var result = [Vector3]()
        for i in 0..<numberOfVertices {
            if vertexVisited[i] {
                result.append(data[i])
            }
        }

        return result
    }

    /// Finds the convex hull and creates the TFace objects.
    private func getConvexFaces() -> [ConvexFace] {
        let faces = convexFaces
        let cellCount = faces.count
        var cells = [ConvexFace]()

        for i in 0..<cellCount {
            let face = facePool[faces[i]];
            var vertices = [Vector3]()
            for j in 0..<3 {
                vertices.append(self.vertices[face.vertices[j]])
            }
            let conFace = ConvexFace(vertices: vertices, normal: isLifted ? Vector3(x:0,y:0,z:0) : face.normal)
            cells.append(conFace)
            face.tag = i;
        }

        for i in 0..<cellCount {
            let face = facePool[faces[i]];

            // Fix the vertex orientation.
            if face.isNormalFlipped {
                var cell = cells[i];

                let tempVert = cells[i].vertices[0];
                cell.vertices[0] = cell.vertices[3 - 1];
                cell.vertices[3 - 1] = tempVert;
                cells[i] = cell

            }
        }

        return cells;
    }


//    /// For 2D only: Returns the result in counter-clockwise order starting with the element with the lowest X value.
//    /// If there are multiple vertices with the same minimum X, then the one with the lowest Y is chosen.
//    private func results2DInOrder<V, TFace: ConvexFace>(data: [V]) -> ConvexHull<TVertex, TFace> where F.TVertex == V{
//        let faces = getConvexFaces();
//    var numPoints = faces.Length;
//    var orderDictionary = new Dictionary<TVertex, TFace>();
//    foreach (var face in faces)
//    orderDictionary.Add(face.Vertices[1], face);
//    var firstPoint = faces[0].Vertices[1];
//    var nextPoint = faces[0].Vertices[0];
//    var orderedPointList = new List<TVertex>();
//    orderedPointList.Add(firstPoint);
//    var orderedFaceList = new List<TFace>();
//    orderedFaceList.Add(faces[1]);
//    var lowestXMinIndex = 0;
//    var k = 0;
//    while (!nextPoint.Equals(firstPoint))
//    {
//    orderedPointList.Add(nextPoint);
//    var nextFace = orderDictionary[nextPoint];
//    orderedFaceList.Add(nextFace);
//    if (nextPoint.Position[0] < orderedPointList[lowestXMinIndex].Position[0]
//    || (nextPoint.Position[0] == orderedPointList[lowestXMinIndex].Position[0]
//    && nextPoint.Position[1] <= orderedPointList[lowestXMinIndex].Position[1]))
//    lowestXMinIndex = k;
//    k++;
//    nextPoint = nextFace.Vertices[0];
//    }
//    TVertex[] points = new TVertex[numPoints];
//    for (int i = 0; i < numPoints; i++)
//    {
//    var j = (i + lowestXMinIndex) % numPoints;
//    points[i] = orderedPointList[j];
//    faces[i] = orderedFaceList[j];
//    }
//    return new ConvexHull<TVertex, TFace>
//    {
//    Points = points,
//    Faces = faces
//    };
//    }


}
