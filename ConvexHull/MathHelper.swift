//
//  MathHelper.swift
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

/// A helper class mostly for normal computation. If convex hulls are computed
/// in higher dimensions, it might be a good idea to add a specific
/// FindNormalVectorND function.
internal class MathHelper {

    /// The matrix pivots
    private var matrixPivots: [Int]

    /// The n d matrix
    private var nDMatrix: [Double]

    /// The n d normal helper vector
    private var nDNormalHelperVector = Vector3(x: 0, y: 0, z: 0)

    /// The nt x
    private var ntX = Vector3(x: 0, y: 0, z: 0)

    /// The nt y
    private var ntY = Vector3(x: 0, y: 0, z: 0)

    /// The nt z
    private var ntZ = Vector3(x: 0, y: 0, z: 0)

    /// The position data
    private var positionData: SimpleList<Double>


    internal init( positions: SimpleList<Double>) {
        self.positionData = positions

        nDMatrix = [Double](repeating:0, count: 9)
        matrixPivots = [0,0,0]
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
        let fi = vertices[0] * 3;
        for i in 0..<3 {
            let n = face.normal[i]
            offset += n * positionData[fi + i];
            centerDistance += n * center[i];
        }
        face.offset = -offset;
        centerDistance -= offset;

        if centerDistance > 0 {
            for i in 0..<3 {
                face.normal[i] = -face.normal[i]
            }
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
        let x = v * 3;
        var distance = f.offset;
        for i in 0..<3 {
            distance += normal[i] * positionData[x + i]
        }
        return distance;
    }

    /// Returns the vector the between vertices.
    internal func vectorBetweenVertices( toIndex: Int, fromIndex: Int) -> Vector3 {
        let u = toIndex * 3, v = fromIndex * 3;
        return Vector3(x: positionData[u] - positionData[v],
                       y: positionData[u + 1] - positionData[v + 1],
                       z: positionData[u + 2] - positionData[v + 2])

    }

    var random = Xoroshiro(seed: (UInt64(arc4random()), UInt64(arc4random())))

    internal func randomOffsetToLift( index: Int, maxHeight: Double) {
        let liftIndex = (index * 3) + 3 - 1;
        let next = random.randomHalfOpen()
        positionData[liftIndex] += 0.0001 * maxHeight * ( next - 0.5)
    }

    /// Finds normal vector of a hyper-plane given by vertices.
    /// Stores the results to normalData.
    private func findNormalVector(vertices: [Int]) -> Vector3 {
        let ntX = vectorBetweenVertices(toIndex: vertices[1], fromIndex: vertices[0]);
        let ntY = vectorBetweenVertices(toIndex: vertices[2], fromIndex: vertices[1]);

        let nx = ntX[1] * ntY[2] - ntX[2] * ntY[1];
        let ny = ntX[2] * ntY[0] - ntX[0] * ntY[2];
        let nz = ntX[0] * ntY[1] - ntX[1] * ntY[0];

        let norm = sqrt(nx * nx + ny * ny + nz * nz);

        let f = 1.0 / norm;
        return Vector3(x: f * nx, y: f * ny, z: f * nz)
    }



    /// Gets the simplex volume. Prior to having enough edge vectors, the method pads the remaining with all
    /// "other numbers". So, yes, this method is not really finding the volume. But a relative volume-like measure. It
    /// uses the magnitude of the determinant as the volume stand-in following the Cayley-Menger theorem.
    internal func getSimplexVolume(edgeVectors: [Vector3], lastIndex: Int, bigNumber: Double) -> Double {
        var A = [Double](repeating: 0, count: 9)
        var index = 0;
        for i in 0..<3 {
            for j in 0..<3 {
                if i <= lastIndex {
                    A[index] = edgeVectors[i][j];
                    index += 1
                } else {
                    let multiplier = index % 2 == 0 ? 1.0 : -1.0

                    A[index] = ( multiplier * Double(index)) / bigNumber
                    index += 1
                }
            }
        }
        // this last term is used for all the vertices in the comparison for the yet determined vertices
        // the idea is to come up with sets of numbers that are orthogonal so that an non-zero value will result
        // and to choose smallish numbers since the choice of vectors will affect what the end volume is.
        // A better way (todo?) is to solve a smaller matrix. However, cases were found in which the obvious smaller vector
        // (the upper left) had too many zeros. So, one would need to find the right subset. Indeed choosing a subset
        // biases the first dimensions of the others. Perhaps a larger volume would be created from a different vertex
        // if another subset of dimensions were used.
        return abs(A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7]
            - A[0] * A[5] * A[7] - A[1] * A[3] * A[8] - A[2] * A[4] * A[6]);
    }


}
