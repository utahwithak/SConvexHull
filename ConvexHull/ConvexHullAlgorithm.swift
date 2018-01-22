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
    private var positions: [Vector3]

    let numberOfVertices: Int

    /// The vertex marks
    private var vertexVisited: [Bool]

    /*
     * The triangulation faces are represented in a single pool for objects that are being reused.
     * This allows for represent the faces as integers and significantly speeds up many computations.
     * - AffectedFaceFlags are used to mark affected faces/
     */
    /// The face pool
    private var facePool = [ConvexFaceInternal]()

    /// The affected face flags
    internal var affectedFaceFlags: [Bool]

    /// Used to track the size of the current hull in the Update/RollbackCenter functions.
    private var convexHullSize = 0

    /// A list of faces that that are not a part of the final convex hull and still need to be processed.
    private var unprocessedFaces = [ConvexFaceInternal]()

    /// A list of faces that form the convex hull.
    private var convexFaces = [Int]()

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
    private var traverseStack = [Int]()


    /// Used for VerticesBeyond for faces that are on the convex hull.
    private let emptyBuffer = [Int]()


    /// Used to determine which vertices are "above" (or "beyond") a face
    private var beyondBuffer = [Int]()


    /// Stores faces that are visible from the current vertex.
    private var affectedFaceBuffer = [Int]()


    /// Stores faces that form a "cone" created by adding new vertex.
    private var coneFaceBuffer = [DeferredFace]()

    /// Stores a list of "singular" (or "generate", "planar", etc.) vertices that cannot be part of the hull.
    private var singularVertices = Set<Int>()


    /// The connector table helps to determine the adjacency of convex faces.
    /// Hashing is used instead of pairwise comparison. This significantly speeds up the computations,
    /// especially for higher dimensions.
    private var connectors = [FaceConnector]()

    /// Helper class for handling math related stuff.
    private var boundingBoxPoints: [[Int]]
    private var indexOfDimensionWithLeastExtremes = 0
    private var minima = [0.0,0.0,0.0]
    private var maxima = [0.0,0.0,0.0]

    // Math Helper
    private var freeFaceIndices = [Int]()
    /// The matrix pivots
    private var matrixPivots = [0,0,0]

    /// The n d matrix
    private var nDMatrix = [Double](repeating:0, count: 9)

    /// The n d normal helper vector
    private var nDNormalHelperVector = Vector3(x: 0, y: 0, z: 0)

    /// The nt x
    private var ntX = Vector3(x: 0, y: 0, z: 0)

    /// The nt y
    private var ntY = Vector3(x: 0, y: 0, z: 0)

    /// The nt z
    private var ntZ = Vector3(x: 0, y: 0, z: 0)


    /// Return the face to the pool for later use.
    internal func depositFace(at faceIndex: Int) {
        let face = facePool[faceIndex]
        face.reset()
        freeFaceIndices.append(faceIndex)
    }

    /// Create a new face and put it in the pool.
    private func createFace() -> ConvexFaceInternal {
        let index = facePool.count
        let face = ConvexFaceInternal(index: index)
        facePool.append(face)
        return face
    }

    public func getFace() -> ConvexFaceInternal {
        if let freeIndex = freeFaceIndices.popLast(){
            return facePool[freeIndex]
        }
        return createFace()
    }

    public static func getConvexHull(with data:[Vector3], planeDistanceTolerance tolerance: Double) -> ConvexHull {

        let ch = ConvexHullAlgorithm(vertices: data, planeDistanceTolerance: tolerance)
        ch.generateConvexHull()
        return ConvexHull(points: ch.hullVertices(data: data), faces: ch.getConvexFaces())
    }

    private init( vertices: [Vector3], planeDistanceTolerance: Double) {
        self.vertices = vertices
        numberOfVertices = vertices.count
        boundingBoxPoints = [[Int]](repeating: [], count: 3)
        vertexVisited = [Bool](repeating: false, count: vertices.count)
        positions = vertices
        affectedFaceFlags = [Bool](repeating: false, count: (3 + 1) * 10)
        updateBuffer = [Int](repeating: 0, count: 3)
        updateIndices = [Int](repeating: 0, count: 3)

        self.planeDistanceTolerance = planeDistanceTolerance
    }

    deinit {
        
    }

    /// Gets/calculates the convex hull. This is
    private func generateConvexHull() {

        // next the bounding box extremes are found. This is used to shift, scale and find the starting simplex.
        findBoundingBoxPoints()
        // the positions are shifted to avoid divide by zero problems
        // and if Delaunay or Voronoi, then the parabola terms are scaled back to match the size of the other coords
        shiftAndScalePositions()
        // Find the (dimension+1) initial points and create the simplexes.
        createInitialSimplex()

        // Now, the main loop. These initial faces of a simplex are replaced and expanded
        // outwards to make the convex hull and faces.
        while let currentFace = unprocessedFaces.first {

            currentVertex = currentFace.furthestVertex

            updateCenter()

            // The affected faces get tagged
            tagAffectedFaces(from: currentFace)

            // Create the cone from the currentVertex and the affected faces horizon.
            if !singularVertices.contains(currentVertex) && createCone() {
                commitCone()
            } else {
                handleSingular()
            }

            // Need to reset the tags
            for bufIndex in affectedFaceBuffer {
                affectedFaceFlags[bufIndex] = false
            }
        }
    }

    /// Finds the bounding box points.
    private func findBoundingBoxPoints() {
        indexOfDimensionWithLeastExtremes = -1
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
                    min = v
                    minIndices.removeAll(keepingCapacity: true)
                    minIndices.append(j)
                } else if (difference > 0) {
                    // you found a solution slightly better than before, clear out those that are no longer on the list and store new value
                    min = v
                    minIndices.removeWhere { min - getCoordinate(vIndex: $0, dimension: i) > planeDistanceTolerance }
                    minIndices.append(j)
                } else if (difference > -planeDistanceTolerance) {
                    //same or almost as good as current limit, so store it
                    minIndices.append(j)
                }

                difference = v - max

                if difference >= planeDistanceTolerance {
                    // you found a better solution than before, clear out the list and store new value
                    max = v
                    maxIndices.removeAll(keepingCapacity: true)
                    maxIndices.append(j)
                } else if difference > 0 {
                    // you found a solution slightly better than before, clear out those that are no longer on the list and store new value
                    max = v
                    maxIndices.removeWhere { min - getCoordinate(vIndex: $0, dimension: i) > planeDistanceTolerance }
                    maxIndices.append(j)
                } else if difference > -planeDistanceTolerance {
                    //same or almost as good as current limit, so store it
                    maxIndices.append(j)
                }
            }
            minima[i] = min
            maxima[i] = max
            minIndices.append(contentsOf: maxIndices)
            if minIndices.count < minNumExtremes {
                minNumExtremes = minIndices.count
                indexOfDimensionWithLeastExtremes = i
            }
            boundingBoxPoints[i] = minIndices
        }
    }

    /// Get a vertex coordinate. In order to reduce speed, all vertex coordinates
    /// have been placed in a single array.
    private func getCoordinate(vIndex: Int, dimension: Int) -> Double {
        let vert = positions[vIndex]
        return dimension == 0 ? vert.x : (dimension == 1 ? vert.y : vert.z)
    }



    /// Shifts and scales the Positions to avoid future errors. This does not alter the original data.
    private func shiftAndScalePositions() {

        var shiftAmount = Vector3.zero

        shiftAmount.x = maxima[0] == minima[0] ? 0 : (maxima[0] - minima[0]) - minima[0]
        shiftAmount.y = maxima[1] == minima[1] ? 0 : (maxima[1] - minima[1]) - minima[1]
        shiftAmount.z = maxima[2] == minima[2] ? 0 : (maxima[2] - minima[2]) - minima[2]

        for i in 0..<positions.count {
            positions[i] += shiftAmount
        }

    }

    /// Find the (dimension+1) initial points and create the simplexes.
    /// Creates the initial simplex of n+1 vertices by using points from the bounding box.
    /// Special care is taken to ensure that the vertices chosen do not result in a degenerate shape
    /// where vertices are collinear (co-planar, etc). This would technically be resolved when additional
    /// vertices are checked in the main loop, but: 1) a degenerate simplex would not eliminate any other
    /// vertices (thus no savings there), 2) the creation of the face normal is prone to error.
    private func createInitialSimplex() {
        var initialPoints = findInitialPoints()
        //create the first faces from (dimension + 1) vertices.
        var faces = [Int](repeating:0, count: 3 + 1)

        for i in 0..<(3 + 1) {
            var vertices = [0,0,0]
            var k = 0
            for j in 0...3 {
                if i != j {
                    vertices[k] = initialPoints[j]
                    k += 1
                }
            }
            let newFace = getFace()

            newFace.vertices = vertices.sorted()

            _ = calculateFacePlane(face: newFace, center: center)
            faces[i] = newFace.index
        }
        // update the adjacency (check all pairs of faces)
        for i in 0..<3 {
            for j in (i + 1)..<(3 + 1) {
                updateAdjacency(l: facePool[faces[i]], r: facePool[faces[j]])
            }
        }


        // Init the vertex beyond buffers.

        for faceIndex in faces  {
            let face = facePool[faceIndex]
            findBeyondVertices(of: face)
            if face.verticesBeyond.count == 0 {
                convexFaces.append(face.index) // The face is on the hull
            } else {
                unprocessedFaces.append(face)
            }
        }


        // Set all vertices to false (unvisited).
        for vertex in initialPoints {
            vertexVisited[vertex] = false
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
        vertexVisited[vertex2] = true
        currentVertex = vertex1
        updateCenter()
        currentVertex = vertex2
        updateCenter()
        var edgeVectors = [Vector3](repeating: .zero, count:3)
        edgeVectors[0] = vectorBetweenVertices(toIndex: vertex2, fromIndex: vertex1)
        // now the remaining vertices are just combined in one big list
        var extremes = boundingBoxPoints.flatMap({ $0})
        // otherwise find the remaining points by maximizing the initial simplex volume
        var index = 1
        while index < 3 && !extremes.isEmpty {
            var bestVertex = -1
            var bestEdgeVector = Vector3.zero
            var maxVolume = Constants.defaultPlaneDistanceTolerance
            for i in stride(from: extremes.count - 1, through: 0, by: -1) {
                // count backwards in order to remove potential duplicates
                let vIndex = extremes[i]
                if initialPoints.contains(vIndex){
                    extremes.remove(at: i)
                } else {
                    edgeVectors[index] = vectorBetweenVertices(toIndex: vIndex, fromIndex: vertex1)
                    let volume = getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber)
                    if maxVolume < volume {
                        maxVolume = volume
                        bestVertex = vIndex
                        bestEdgeVector = edgeVectors[index]
                    }
                }
            }
            if let index = extremes.index(of: bestVertex) {
                extremes.remove(at: index)
            }
            if bestVertex == -1 {
                break
            }
            initialPoints.append(bestVertex)
            edgeVectors[index] = bestEdgeVector
            index += 1
            currentVertex = bestVertex
            updateCenter()
        }
        // hmm, there are not enough points on the bounding box to make a simplex. It is rare but entirely possibly.
        // As an extreme, the bounding box can be made in n dimensions from only 2 unique points. When we can't find
        // enough unique points, we start again with ALL the vertices. The following is a near replica of the code
        // above, but instead of extremes, we consider "allVertices".
        if initialPoints.count <= 3 {
            var allVertices = [Int](0..<numberOfVertices)
            while index < 3 && !allVertices.isEmpty {
                var bestVertex = -1
                var bestEdgeVector = Vector3.zero
                var maxVolume = 0.0
                for i in stride(from: allVertices.count - 1, through: 0, by: -1) {
                    // count backwards in order to remove potential duplicates
                    let vIndex = allVertices[i]
                    if initialPoints.contains(vIndex) {
                        allVertices.remove(at: i)
                    } else {
                        edgeVectors[index] = vectorBetweenVertices(toIndex: vIndex, fromIndex: vertex1)
                        let volume = getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber)
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
                    break
                }
                initialPoints.append(bestVertex)
                edgeVectors[index] = bestEdgeVector
                index += 1
                currentVertex = bestVertex
                updateCenter()
            }
        }

        return initialPoints
    }


    /// Check if 2 faces are adjacent and if so, update their AdjacentFaces array.
    private func updateAdjacency(l: ConvexFaceInternal, r: ConvexFaceInternal) {
        let lv = l.vertices
        let rv = r.vertices

        // reset marks on the 1st face
        lv.forEach { vertexVisited[$0] = false}

        // mark all vertices on the 2nd face
        rv.forEach({ vertexVisited[$0] = true })

        var i = 0

        // find the 1st false index
        for k in 0..<lv.count {
            if !vertexVisited[lv[k]] {
                i = k
                break
            }
        }
        // no vertex was marked
        if i == 3 {
            return
        }


        // check if only 1 vertex wasn't marked
        for j in (i + 1)..<lv.count {
            if !vertexVisited[lv[j]] {
                return
            }
        }

        // if we are here, the two faces share an edge
        l.adjacentFaces[i] = r.index

        // update the adj. face on the other face - find the vertex that remains marked
        lv.forEach { vertexVisited[$0] = false}
        for k in 0..<rv.count {
            if vertexVisited[rv[k]] {
                i = k
                break
            }
        }
        r.adjacentFaces[i] = l.index
    }


    /// Used in the "initialization" code.
    private func findBeyondVertices(of face: ConvexFaceInternal) {
        maxDistance = -Double.greatestFiniteMagnitude
        furthestVertex = 0
        for i in 0..<numberOfVertices where !vertexVisited[i] {
            isBeyond(face: face, beyondVertices: &face.verticesBeyond, v: i)
        }

        face.furthestVertex = furthestVertex
    }


    /// Check whether the vertex v is beyond the given face. If so, add it to beyondVertices.
    private func isBeyond(face: ConvexFaceInternal, beyondVertices: inout [Int], v: Int) {
        let distance = getVertexDistance(v: v, f: face)
        if distance >= planeDistanceTolerance {
            if distance > maxDistance {
                // If it's within the tolerance distance, use the lex. larger point
                if distance - maxDistance < planeDistanceTolerance { // todo: why is this LexCompare necessary. Would seem to favor x over y over z (etc.)?
                    if lexCompare(u: v, v: furthestVertex) > 0 {
                        maxDistance = distance
                        furthestVertex = v
                    }
                } else {
                    maxDistance = distance
                    furthestVertex = v
                }
            }
            beyondVertices.append(v)
        }
    }

    /// Compares the values of two vertices. The return value (-1, 0 or +1) are found
    /// by first checking the first coordinate and then progressing through the rest.
    /// In this way {2, 8} will be a "-1" (less than) {3, 1}.
    private func lexCompare(u: Int, v: Int) -> Int {

        let a = positions[u]
        let b = positions[v]
        if a.x < b.x {
            return -1
        } else if a.x > b.x {
            return 1
        } else if a.y < b.y {
            return -1
        } else if a.y > b.y {
            return 1
        } else if a.z < b.z {
            return -1
        } else if a.z > b.z {
            return 1
        }

        return 0
    }

    /// Recalculates the centroid of the current hull.
    private func updateCenter() {

        center *= Double(convexHullSize)

        convexHullSize += 1
        let f = 1.0 / Double(convexHullSize)
        let co = currentVertex
        center = f * (center + positions[co])

    }


    /// Tags all faces seen from the current vertex with 1.
    private func tagAffectedFaces(from currentFace: ConvexFaceInternal) {
        affectedFaceBuffer.removeAll(keepingCapacity: true)
        affectedFaceBuffer.append(currentFace.index)
        traverseAffectedFaces(from: currentFace.index)
    }


    /// Recursively traverse all the relevant faces.
    private func traverseAffectedFaces(from currentFace: Int) {
        traverseStack.removeAll(keepingCapacity: true)
        traverseStack.append(currentFace)
        affectedFaceFlags[currentFace] = true
        while let toVisit = traverseStack.popLast() {
            let top = facePool[toVisit]
            for i in 0..<3 {
                let adjFace = top.adjacentFaces[i]

                if !affectedFaceFlags[adjFace] && getVertexDistance(v: currentVertex, f: facePool[adjFace]) >= planeDistanceTolerance {
                    affectedFaceBuffer.append(adjFace)
                    affectedFaceFlags[adjFace] = true
                    traverseStack.append(adjFace)
                }
            }
        }
    }

    /// Connect faces using a connector.
    private func connectFace(with connector: FaceConnector) {
        let hash = connector.hashCode

        for (index, current) in connectors.enumerated() where current.hashCode == hash {
            if connector.hashCode == current.hashCode && connector.v0 == current.v0 && connector.v1 == current.v1 {
                connectors.remove(at: index)
                current.face.adjacentFaces[current.edgeIndex] = connector.face.index;
                connector.face.adjacentFaces[connector.edgeIndex] = current.face.index;
                return
            }
        }
        connectors.append(connector)
    }

    /// Removes the faces "covered" by the current vertex and adds the newly created ones.
    private func createCone() -> Bool {
        let currentVertexIndex = currentVertex
        coneFaceBuffer.removeAll(keepingCapacity: true)

        for fIndex in 0..<affectedFaceBuffer.count {
            let oldFaceIndex = affectedFaceBuffer[fIndex]
            let oldFace = facePool[oldFaceIndex]

            // Find the faces that need to be updated
            var updateCount = 0
            for i in 0..<3 {
                let af = oldFace.adjacentFaces[i]
                if !affectedFaceFlags[af] {
                    updateBuffer[updateCount] = af
                    updateIndices[updateCount] = i
                    updateCount += 1
                }
            }

            for i in 0..<updateCount {
                let adjacentFace = facePool[updateBuffer[i]]

                var oldFaceAdjacentIndex = 0
                let adjFaceAdjacency = adjacentFace.adjacentFaces
                for j in 0..<adjFaceAdjacency.count {
                    if (oldFaceIndex == adjFaceAdjacency[j]) {
                        oldFaceAdjacentIndex = j
                        break
                    }
                }

                let forbidden = updateIndices[i] // Index of the face that corresponds to this adjacent face

                var oldVertexIndex = 0

                let newFace = getFace()

                for j in 0..<3 {
                    newFace.vertices[j] = oldFace.vertices[j]
                }
                oldVertexIndex = newFace.vertices[forbidden]

                var orderedPivotIndex = 0

                // correct the ordering
                if currentVertexIndex < oldVertexIndex {
                    orderedPivotIndex = 0
                    for j in stride(from: forbidden - 1,through: 0, by: -1) {
                        if newFace.vertices[j] > currentVertexIndex {
                            newFace.vertices[j + 1] = newFace.vertices[j]
                        } else {
                            orderedPivotIndex = j + 1
                            break
                        }
                    }
                } else {
                    orderedPivotIndex = 3 - 1
                    for j in (forbidden + 1)..<3 {
                        if newFace.vertices[j] < currentVertexIndex {
                            newFace.vertices[j - 1] = newFace.vertices[j]
                        } else {
                            orderedPivotIndex = j - 1
                            break
                        }
                    }
                }

                newFace.vertices[orderedPivotIndex] = currentVertex

                if !calculateFacePlane(face: newFace, center: center) {
                    return false
                }
                let deferredFace = DeferredFace(face: newFace, pivot: adjacentFace, oldFace: oldFace, faceIndex: orderedPivotIndex, pivotIndex: oldFaceAdjacentIndex)

                coneFaceBuffer.append(deferredFace)

            }
        }

        return true
    }

    /// Commits a cone and adds a vertex to the convex hull.
    private func commitCone() {
        // Fill the adjacency.
        for face in coneFaceBuffer {

            let newFace = face.face
            let adjacentFace = face.pivot
            let oldFace = face.oldFace

            let orderedPivotIndex = face.faceIndex

            newFace.adjacentFaces[orderedPivotIndex] = adjacentFace.index
            adjacentFace.adjacentFaces[face.pivotIndex] = newFace.index

            // let there be a connection.
            for j in 0..<3 {
                if (j == orderedPivotIndex) {
                    continue
                }

                let connector = FaceConnector(face: newFace, edgeIndex: j)
                connectFace(with: connector)
            }

            // the id adjacent face on the hull? If so, we can use simple method to find beyond vertices.
            if adjacentFace.verticesBeyond.count == 0 {
                findBeyondVertices(face: newFace, beyond: oldFace.verticesBeyond)
            } else if (adjacentFace.verticesBeyond.count < oldFace.verticesBeyond.count) { // it is slightly more effective if the face with the lower number of beyond vertices comes first.
                findBeyondVertices(face: newFace, beyond: adjacentFace.verticesBeyond, beyond1: oldFace.verticesBeyond)
            } else {
                findBeyondVertices(face: newFace, beyond: oldFace.verticesBeyond, beyond1: adjacentFace.verticesBeyond)
            }

            // This face will definitely lie on the hull
            if newFace.verticesBeyond.isEmpty {
                convexFaces.append(newFace.index)
                if let index = unprocessedFaces.index(where: { $0 === newFace }) {
                    unprocessedFaces.remove(at: index)
                }
                newFace.verticesBeyond.removeAll(keepingCapacity: true)
            } else {
                unprocessedFaces.append(newFace)
            }

        }

        // Recycle the affected faces.
        for fIndex in 0..<affectedFaceBuffer.count {
            let faceIndex = affectedFaceBuffer[fIndex]
            if let index = unprocessedFaces.index(where: { $0 === facePool[faceIndex]}) {
                unprocessedFaces.remove(at: index)
            }
            depositFace(at: faceIndex)
        }
    }

    /// Handles singular vertex.
    private func handleSingular() {
        rollbackCenter()
        singularVertices.insert(currentVertex)

        // This means that all the affected faces must be on the hull and that all their "vertices beyond" are singular.
        for fIndex in 0..<affectedFaceBuffer.count {
            let face = facePool[affectedFaceBuffer[fIndex]]
            let vb = face.verticesBeyond
            for i in 0..<vb.count {
                singularVertices.insert(vb[i])
            }

            convexFaces.append(face.index)
            if let index = unprocessedFaces.index(where: { $0 === face} ) {
                unprocessedFaces.remove(at: index)
            }
            face.verticesBeyond.removeAll(keepingCapacity: true)
        }
    }

    /// Removes the last vertex from the center.
    private func rollbackCenter() {
        center *= Double(convexHullSize)

        convexHullSize -= 1
        let f = convexHullSize > 0 ? 1.0 / Double(convexHullSize) : 0.0
        let co = currentVertex

        center = f * (center - positions[co])

    }

    /// Used by update faces.
    private func findBeyondVertices(face: ConvexFaceInternal, beyond: [Int] , beyond1: [Int] ) {
        maxDistance = Double.greatestFiniteMagnitude * -1
        furthestVertex = 0
        var v = 0

        for i in 0..<beyond1.count {
            vertexVisited[beyond1[i]] = true
        }

        vertexVisited[currentVertex] = false
        for i in 0..<beyond.count {
            v = beyond[i]
            if v == currentVertex{
                continue
            }
            vertexVisited[v] = false
            isBeyond(face: face, beyondVertices: &beyondBuffer, v: v)
        }

        for i in 0..<beyond1.count {
            v = beyond1[i]
            if vertexVisited[v] {
                isBeyond(face: face, beyondVertices: &beyondBuffer, v: v)
            }
        }

        face.furthestVertex = furthestVertex

        // Pull the old switch a roo (switch the face beyond buffers)
        var temp = face.verticesBeyond
        face.verticesBeyond = beyondBuffer
        if (temp.count > 0){
            temp.removeAll(keepingCapacity: true)
        }
        beyondBuffer = temp
    }

    /// Finds the beyond vertices.
    private func findBeyondVertices(face: ConvexFaceInternal , beyond: [Int]){

        maxDistance = -Double.greatestFiniteMagnitude
        furthestVertex = 0
        var v = 0

        for i in 0..<beyond.count {
            v = beyond[i]
            if v == currentVertex {
                continue
            }
            isBeyond(face: face, beyondVertices: &beyondBuffer, v: v)
        }

        face.furthestVertex = furthestVertex

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
            let face = facePool[faces[i]]
            var vertices = [Vector3]()
            for j in 0..<3 {
                vertices.append(self.vertices[face.vertices[j]])
            }
            let conFace = ConvexFace(vertices: vertices, normal: face.normal)
            cells.append(conFace)
            face.tag = i
        }

        for i in 0..<cellCount {
            let face = facePool[faces[i]]

            // Fix the vertex orientation.
            if face.isNormalFlipped {
                var cell = cells[i]

                let tempVert = cells[i].vertices[0]
                cell.vertices[0] = cell.vertices[3 - 1]
                cell.vertices[3 - 1] = tempVert
                cells[i] = cell

            }
        }

        return cells
    }

    /// Calculates the normal and offset of the hyper-plane given by the face's vertices.
    internal func calculateFacePlane(face: ConvexFaceInternal, center: Vector3) -> Bool {
        var vertices = face.vertices;

        face.normal = findNormalVector(vertices: vertices);

        if face.normal.x.isNaN {
            return false
        }

        var offset = 0.0;
        var centerDistance = 0.0;
        let fi = vertices[0]
        let faceNorm = face.normal
        offset += faceNorm.x * positions[fi].x
        centerDistance += faceNorm.x * center.x
        offset += faceNorm.y * positions[fi].y
        centerDistance += faceNorm.y * center.y
        offset += faceNorm.z * positions[fi].z
        centerDistance += faceNorm.z * center.z

        face.offset = -offset;
        centerDistance -= offset;

        if centerDistance > 0 {

            face.normal = face.normal * -1
            face.offset = offset;
            face.isNormalFlipped = true
        } else {
            face.isNormalFlipped = false
        }

        return true;
    }

    /// Check if the vertex is "visible" from the face.
    /// The vertex is "over face" if the return value is > Constants.PlaneDistanceTolerance.
    internal func getVertexDistance(v: Int, f: ConvexFaceInternal) -> Double {
        let normal = f.normal

        var distance = f.offset;

        distance += normal.x * positions[v].x
        distance += normal.y * positions[v].y
        distance += normal.z * positions[v].z

        return distance;
    }

    /// Returns the vector the between vertices.
    internal func vectorBetweenVertices( toIndex: Int, fromIndex: Int) -> Vector3 {

        return positions[toIndex] - positions[fromIndex]

    }

    var random = Xoroshiro(seed: (UInt64(arc4random()), UInt64(arc4random())))

    internal func randomOffsetToLift( index: Int, maxHeight: Double) {
        let next = random.randomHalfOpen()
        positions[index].y += 0.0001 * maxHeight * ( next - 0.5)
    }

    /// Finds normal vector of a hyper-plane given by vertices.
    /// Stores the results to normalData.
    private func findNormalVector(vertices: [Int]) -> Vector3 {
        let ntX = vectorBetweenVertices(toIndex: vertices[1], fromIndex: vertices[0]);
        let ntY = vectorBetweenVertices(toIndex: vertices[2], fromIndex: vertices[1]);

        let nx = ntX.y * ntY.z - ntX.z * ntY.y
        let ny = ntX.z * ntY.x - ntX.x * ntY.z
        let nz = ntX.x * ntY.y - ntX.y * ntY.x

        let norm = sqrt(nx * nx + ny * ny + nz * nz);

        let f = 1.0 / norm;
        return Vector3(x: f * nx, y: f * ny, z: f * nz)
    }

    /// Gets the simplex volume. Prior to having enough edge vectors, the method pads the remaining with all
    /// "other numbers". So, yes, this method is not really finding the volume. But a relative volume-like measure. It
    /// uses the magnitude of the determinant as the volume stand-in following the Cayley-Menger theorem.
    internal func getSimplexVolume(edgeVectors: [Vector3], lastIndex: Int, bigNumber: Double) -> Double {

        let a0 = 0 <= lastIndex ? edgeVectors[0].x :  0
        let a1 = 0 <= lastIndex ? edgeVectors[0].y : -1 / bigNumber
        let a2 = 0 <= lastIndex ? edgeVectors[0].z :  2 / bigNumber
        let a3 = 1 <= lastIndex ? edgeVectors[1].x : -3 / bigNumber
        let a4 = 1 <= lastIndex ? edgeVectors[1].y :  4 / bigNumber
        let a5 = 1 <= lastIndex ? edgeVectors[1].z : -5 / bigNumber
        let a6 = 2 <= lastIndex ? edgeVectors[2].x :  6 / bigNumber
        let a7 = 2 <= lastIndex ? edgeVectors[2].y : -7 / bigNumber
        let a8 = 2 <= lastIndex ? edgeVectors[2].z :  8 / bigNumber


        // this last term is used for all the vertices in the comparison for the yet determined vertices
        // the idea is to come up with sets of numbers that are orthogonal so that an non-zero value will result
        // and to choose smallish numbers since the choice of vectors will affect what the end volume is.
        // A better way (todo?) is to solve a smaller matrix. However, cases were found in which the obvious smaller vector
        // (the upper left) had too many zeros. So, one would need to find the right subset. Indeed choosing a subset
        // biases the first dimensions of the others. Perhaps a larger volume would be created from a different vertex
        // if another subset of dimensions were used.
        return abs(a0 * a4 * a8 + a1 * a5 * a6 + a2 * a3 * a7
            - a0 * a5 * a7 - a1 * a3 * a8 - a2 * a4 * a6);
    }



}
