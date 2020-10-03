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

class ConvexHullAlgorithm {

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
    internal var affectedFaceFlags =  [Bool]()

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
        affectedFaceFlags.append(false)
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

    private init(vertices: [Vector3], planeDistanceTolerance: Double) {
        self.vertices = vertices
        numberOfVertices = vertices.count
        boundingBoxPoints = [[Int]](repeating: [], count: 3)
        vertexVisited = [Bool](repeating: false, count: vertices.count)
        positions = vertices
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
        let initialPoints = findInitialPoints()
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
            vertices.sort()
            newFace.vert0 = vertices[0]
            newFace.vert1 = vertices[1]
            newFace.vert2 = vertices[2]

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
            if face.verticesBeyond.isEmpty {
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
        edgeVectors[0] = positions[vertex2] - positions[vertex1]
        // now the remaining vertices are just combined in one big list
        var extremes = boundingBoxPoints.flatMap({ $0})
        // otherwise find the remaining points by maximizing the initial simplex volume
        var index = 1
        while index < 3 && !extremes.isEmpty {
            var bestVertex = -1
            var bestEdgeVector = Vector3.zero
            var maxVolume = ConvexHull.defaultPlaneDistanceTolerance
            for i in stride(from: extremes.count - 1, through: 0, by: -1) {
                // count backwards in order to remove potential duplicates
                let vIndex = extremes[i]
                if initialPoints.contains(vIndex){
                    extremes.remove(at: i)
                } else {
                    edgeVectors[index] = positions[vIndex] - positions[vertex1]
                    let volume = getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber)
                    if maxVolume < volume {
                        maxVolume = volume
                        bestVertex = vIndex
                        bestEdgeVector = edgeVectors[index]
                    }
                }
            }
            if let index = extremes.firstIndex(of: bestVertex) {
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
                        edgeVectors[index] = positions[vIndex] - positions[vertex1]
                        let volume = getSimplexVolume(edgeVectors: edgeVectors, lastIndex: index, bigNumber: bigNumber)
                        if maxVolume < volume {
                            maxVolume = volume
                            bestVertex = vIndex
                            bestEdgeVector = edgeVectors[index]
                        }
                    }
                }
                if let index = allVertices.firstIndex(of: bestVertex) {
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

        vertexVisited[l.vert0] = false
        vertexVisited[l.vert1] = false
        vertexVisited[l.vert2] = false

        vertexVisited[r.vert0] = true
        vertexVisited[r.vert1] = true
        vertexVisited[r.vert2] = true


        if !vertexVisited[l.vert0] && vertexVisited[l.vert1] && vertexVisited[l.vert2] {
            l.adj0 = r.index
        } else if !vertexVisited[l.vert1] && vertexVisited[l.vert0] && vertexVisited[l.vert2] {
            l.adj1 = r.index
        } else if !vertexVisited[l.vert2] && vertexVisited[l.vert0] && vertexVisited[l.vert1]  {
            l.adj2 = r.index
        } else {
            return
        }

        // update the adj. face on the other face - find the vertex that remains marked
        vertexVisited[l.vert0] = false
        vertexVisited[l.vert1] = false
        vertexVisited[l.vert2] = false

        if vertexVisited[r.vert0] {
            r.adj0 = l.index
        } else if vertexVisited[r.vert1] {
            r.adj1 = l.index
        } else if vertexVisited[r.vert2] {
            r.adj2 = l.index
        }

    }


    /// Used in the "initialization" code.
    private func findBeyondVertices(of face: ConvexFaceInternal) {
        maxDistance = -Double.greatestFiniteMagnitude
        furthestVertex = 0
        for i in 0..<numberOfVertices where !vertexVisited[i] {
            if isBeyond(face: face, v: i) {
                face.verticesBeyond.append(i)
            }
        }

        face.furthestVertex = furthestVertex
    }


    /// Check whether the vertex v is beyond the given face. If so, add it to beyondVertices.
    private func isBeyond(face: ConvexFaceInternal, v: Int) -> Bool {
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
            return true
        }
        return false

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

            if !affectedFaceFlags[top.adj0] && getVertexDistance(v: currentVertex, f: facePool[top.adj0]) >= planeDistanceTolerance {
                affectedFaceBuffer.append(top.adj0)
                affectedFaceFlags[top.adj0] = true
                traverseStack.append(top.adj0)
            }

            if !affectedFaceFlags[top.adj1] && getVertexDistance(v: currentVertex, f: facePool[top.adj1]) >= planeDistanceTolerance {
                affectedFaceBuffer.append(top.adj1)
                affectedFaceFlags[top.adj1] = true
                traverseStack.append(top.adj1)
            }

            if !affectedFaceFlags[top.adj2] && getVertexDistance(v: currentVertex, f: facePool[top.adj2]) >= planeDistanceTolerance {
                affectedFaceBuffer.append(top.adj2)
                affectedFaceFlags[top.adj2] = true
                traverseStack.append(top.adj2)
            }

        }
    }

    /// Connect faces using a connector.
    private func connectFace(with connector: FaceConnector) {
        let hash = connector.hashCode
        let v0 = connector.v0
        let v1 = connector.v1
        if let index = connectors.firstIndex(where: { $0.hashCode == hash && v0 == $0.v0 && v1 == $0.v1 }) {
            let current = connectors.remove(at: index)

            if current.edgeIndex == 0 {
                current.face.adj0 = connector.face.index
            } else if current.edgeIndex == 1 {
                current.face.adj1 = connector.face.index
            } else {
                current.face.adj2 = connector.face.index
            }

            if connector.edgeIndex == 0 {
                connector.face.adj0 = current.face.index
            } else if connector.edgeIndex == 1 {
                connector.face.adj1 = current.face.index
            } else {
                connector.face.adj2 = current.face.index
            }

        } else {
            connectors.append(connector)
        }
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

            if !affectedFaceFlags[oldFace.adj0] {
                updateBuffer[updateCount] = oldFace.adj0
                updateIndices[updateCount] = 0
                updateCount += 1
            }

            if !affectedFaceFlags[oldFace.adj1] {
                updateBuffer[updateCount] = oldFace.adj1
                updateIndices[updateCount] = 1
                updateCount += 1
            }

            if !affectedFaceFlags[oldFace.adj2] {
                updateBuffer[updateCount] = oldFace.adj2
                updateIndices[updateCount] = 2
                updateCount += 1
            }


            for i in 0..<updateCount {
                let adjacentFace = facePool[updateBuffer[i]]

                let oldFaceAdjacentIndex: Int

                if oldFaceIndex == adjacentFace.adj0 {
                    oldFaceAdjacentIndex = 0
                } else if oldFaceIndex == adjacentFace.adj1 {
                    oldFaceAdjacentIndex = 1
                } else if oldFaceIndex == adjacentFace.adj2 {
                    oldFaceAdjacentIndex = 2
                } else {
                    oldFaceAdjacentIndex = 0
                }

                let forbidden = updateIndices[i] // Index of the face that corresponds to this adjacent face


                let newFace = getFace()
                newFace.vert0 = oldFace.vert0
                newFace.vert1 = oldFace.vert1
                newFace.vert2 = oldFace.vert2

                let oldVertexIndex: Int
                if forbidden == 0 {
                    oldVertexIndex = newFace.vert0
                } else if forbidden == 1 {
                    oldVertexIndex = newFace.vert1
                } else {
                    oldVertexIndex = newFace.vert2
                }

                var orderedPivotIndex = 0

                // correct the ordering
                if currentVertexIndex < oldVertexIndex {
                    orderedPivotIndex = 0
                    for j in stride(from: forbidden - 1,through: 0, by: -1) {
                        if newFace[j] > currentVertexIndex {
                            newFace[j + 1] = newFace[j]
                        } else {
                            orderedPivotIndex = j + 1
                            break
                        }
                    }
                } else {
                    orderedPivotIndex = 3 - 1
                    for j in (forbidden + 1)..<3 {
                        if newFace[j] < currentVertexIndex {
                            newFace[j - 1] = newFace[j]
                        } else {
                            orderedPivotIndex = j - 1
                            break
                        }
                    }
                }
                if orderedPivotIndex == 0 {
                    newFace.vert0 = currentVertex
                } else if orderedPivotIndex == 1 {
                    newFace.vert1 = currentVertex
                } else {
                    newFace.vert2 = currentVertex
                }

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

            newFace.set(adj: orderedPivotIndex, to: adjacentFace.index)
            adjacentFace.set(adj: face.pivotIndex, to: newFace.index)

            // let there be a connection.
            if orderedPivotIndex == 0 {
                connectFace(with: FaceConnector(face: newFace, edgeIndex: 1))
                connectFace(with: FaceConnector(face: newFace, edgeIndex: 2))
            } else if orderedPivotIndex == 1 {
                connectFace(with: FaceConnector(face: newFace, edgeIndex: 0))
                connectFace(with: FaceConnector(face: newFace, edgeIndex: 2))
            } else {
                connectFace(with: FaceConnector(face: newFace, edgeIndex: 0))
                connectFace(with: FaceConnector(face: newFace, edgeIndex: 1))
            }

            // the id adjacent face on the hull? If so, we can use simple method to find beyond vertices.
            if adjacentFace.verticesBeyond.isEmpty {
                findBeyondVertices(face: newFace, beyond: oldFace.verticesBeyond)
            } else if (adjacentFace.verticesBeyond.count < oldFace.verticesBeyond.count) { // it is slightly more effective if the face with the lower number of beyond vertices comes first.
                findBeyondVertices(face: newFace, beyond: adjacentFace.verticesBeyond, beyond1: oldFace.verticesBeyond)
            } else {
                findBeyondVertices(face: newFace, beyond: oldFace.verticesBeyond, beyond1: adjacentFace.verticesBeyond)
            }

            // This face will definitely lie on the hull
            if newFace.verticesBeyond.isEmpty {
                convexFaces.append(newFace.index)
                if let index = unprocessedFaces.firstIndex(where: { $0 === newFace }) {
                    unprocessedFaces.remove(at: index)
                }
                newFace.verticesBeyond.removeAll(keepingCapacity: true)
            } else {
                unprocessedFaces.append(newFace)
            }

        }

        // Recycle the affected faces.
        for faceIndex in affectedFaceBuffer {

            if let index = unprocessedFaces.firstIndex(where: { $0 === facePool[faceIndex]}) {
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
            if let index = unprocessedFaces.firstIndex(where: { $0 === face} ) {
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

        if !face.verticesBeyond.isEmpty {
            face.verticesBeyond.removeAll(keepingCapacity: true)
        }
        for v in beyond1 {
            vertexVisited[v] = true
        }

        vertexVisited[currentVertex] = false
        for v in beyond where v != currentVertex {

            vertexVisited[v] = false
            if isBeyond(face: face, v: v) {
                face.verticesBeyond.append(v)
            }
        }

        for v in beyond1 where vertexVisited[v] && isBeyond(face: face, v: v) {
            face.verticesBeyond.append(v)
        }

        face.furthestVertex = furthestVertex
    }

    /// Finds the beyond vertices.
    private func findBeyondVertices(face: ConvexFaceInternal , beyond: [Int]) {
        if !face.verticesBeyond.isEmpty {
            face.verticesBeyond.removeAll(keepingCapacity: true)
        }

        maxDistance = -Double.greatestFiniteMagnitude
        furthestVertex = 0

        for v in beyond where v != currentVertex && isBeyond(face: face, v: v) {
            face.verticesBeyond.append(v)
        }

        face.furthestVertex = furthestVertex

    }


    /// Gets the hull vertices.
    private func hullVertices(data: [Vector3]) -> [Vector3] {
        let cellCount = convexFaces.count

        for i in 0..<numberOfVertices {
            vertexVisited[i] = false
        }

        for i in 0..<cellCount {
            let face = facePool[convexFaces[i]]
            vertexVisited[face.vert0] = true
            vertexVisited[face.vert1] = true
            vertexVisited[face.vert2] = true
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
            let conFace = ConvexFace(vertices: [self.vertices[face.vert0], self.vertices[face.vert1], self.vertices[face.vert2]], normal: face.normal)
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
    func calculateFacePlane(face: ConvexFaceInternal, center: Vector3) -> Bool {

        face.normal = findNormalVector(a: face.vert0, b: face.vert1, c: face.vert2)

        if face.normal.x.isNaN {
            return false
        }

        var offset = 0.0;
        var centerDistance = 0.0;
        let fi = face.vert0
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
    func getVertexDistance(v: Int, f: ConvexFaceInternal) -> Double {
        let normal = f.normal

        var distance = f.offset;

        distance += normal.x * positions[v].x
        distance += normal.y * positions[v].y
        distance += normal.z * positions[v].z

        return distance;
    }

    var random = Xoroshiro(seed: (UInt64(arc4random()), UInt64(arc4random())))

    internal func randomOffsetToLift( index: Int, maxHeight: Double) {
        let next = random.randomHalfOpen()
        positions[index].y += 0.0001 * maxHeight * ( next - 0.5)
    }

    /// Finds normal vector of a hyper-plane given by vertices.
    /// Stores the results to normalData.
    private func findNormalVector(a: Int, b: Int, c: Int) -> Vector3 {
        let ntX = positions[b] - positions[a]
        let ntY = positions[c] - positions[b]

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
