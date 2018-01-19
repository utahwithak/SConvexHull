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

    /// The dimension
    private let dimension: Int

    /// The matrix pivots
    private var matrixPivots : [Int]

    /// The n d matrix
    private var nDMatrix: [Double]

    /// The n d normal helper vector
    private var nDNormalHelperVector: [Double]

    /// The nt x
    private var ntX: [Double]

    /// The nt y
    private var ntY: [Double]

    /// The nt z
    private var ntZ: [Double]

    /// The position data
    private var positionData: SimpleList<Double>


    internal init(dimension: Int,  positions: SimpleList<Double>) {
        self.positionData = positions
        self.dimension = dimension

        ntX = [Double](repeating:0, count: dimension)
        ntY = [Double](repeating:0, count: dimension)
        ntZ = [Double](repeating:0, count: dimension)

        nDNormalHelperVector = [Double](repeating:0, count: dimension)
        nDMatrix = [Double](repeating:0, count: dimension * dimension)
        matrixPivots = [Int](repeating:0, count: dimension)
    }

    /// Calculates the normal and offset of the hyper-plane given by the face's vertices.
    internal func calculateFacePlane(face: ConvexFaceInternal, center: [Double]) -> Bool {
        var vertices = face.vertices;

        findNormalVector(vertices: vertices, normalData: &face.normal);

        if face.normal[0].isNaN {
            return false
        }

        var offset = 0.0;
        var centerDistance = 0.0;
        let fi = vertices[0] * dimension;
        for i in 0..<dimension {
            let n = face.normal[i]
            offset += n * positionData[fi + i];
            centerDistance += n * center[i];
        }
        face.offset = -offset;
        centerDistance -= offset;

        if centerDistance > 0 {
            for i in 0..<dimension {
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
    /// The vertex is "over face" if the return value is &gt; Constants.PlaneDistanceTolerance.
    internal func getVertexDistance(v: Int, f: ConvexFaceInternal ) -> Double {
        let normal = f.normal
        let x = v * dimension;
        var distance = f.offset;
        for i in 0..<normal.count {
            distance += normal[i] * positionData[x + i]
        }
        return distance;
    }

    /// Returns the vector the between vertices.
    internal func vectorBetweenVertices(toIndex: Int, fromIndex: Int) -> [Double] {
        var target = [Double](repeating: 0, count: dimension)
        vectorBetweenVertices(toIndex: toIndex, fromIndex: fromIndex, target: &target);
        return target;
    }

    /// Returns the vector the between vertices.
    private func vectorBetweenVertices( toIndex: Int, fromIndex: Int, target: inout [Double]) {
        let u = toIndex * dimension, v = fromIndex * dimension;
        for i in 0..<dimension {
            target[i] = positionData[u + i] - positionData[v + i];
        }
    }

    var random = Xoroshiro(seed: (UInt64(arc4random()), UInt64(arc4random())))

    internal func randomOffsetToLift( index: Int, maxHeight: Double) {
        let liftIndex = (index * dimension) + dimension - 1;
        let next = random.randomHalfOpen()
        positionData[liftIndex] += 0.0001 * maxHeight * ( next - 0.5)
    }

    /// Finds normal vector of a hyper-plane given by vertices.
    /// Stores the results to normalData.
    private func findNormalVector(vertices: [Int], normalData: inout [Double]) {
        switch dimension {
        case 2:
            findNormalVector2D(vertices: vertices, normal: &normalData);
        case 3:
            findNormalVector3D(vertices: vertices, normal: &normalData);
        case 4:
            findNormalVector4D(vertices: vertices, normal: &normalData);
        default:
            findNormalVectorND(vertices: vertices, normal: &normalData);
        }
    }

    /// Finds 2D normal vector.
    private func findNormalVector2D(vertices: [Int], normal: inout [Double]) {

        vectorBetweenVertices(toIndex: vertices[1], fromIndex: vertices[0], target: &ntX);

        let nx = -ntX[1];
        let ny = ntX[0];

        let norm = sqrt(nx * nx + ny * ny);

        let f = 1.0 / norm;
        normal[0] = f * nx;
        normal[1] = f * ny;
    }

    /// Finds 3D normal vector
    private func findNormalVector3D( vertices: [Int], normal: inout [Double]) {
        vectorBetweenVertices(toIndex: vertices[1], fromIndex: vertices[0], target: &ntX);
        vectorBetweenVertices(toIndex: vertices[2], fromIndex: vertices[1], target: &ntY);

        let nx = ntX[1] * ntY[2] - ntX[2] * ntY[1];
        let ny = ntX[2] * ntY[0] - ntX[0] * ntY[2];
        let nz = ntX[0] * ntY[1] - ntX[1] * ntY[0];

        let norm = sqrt(nx * nx + ny * ny + nz * nz);

        let f = 1.0 / norm;
        normal[0] = f * nx;
        normal[1] = f * ny;
        normal[2] = f * nz;
    }

    /// Finds 4D normal vector
    private func findNormalVector4D(vertices: [Int], normal: inout [Double]) {
        vectorBetweenVertices(toIndex: vertices[1], fromIndex: vertices[0], target: &ntX);
        vectorBetweenVertices(toIndex: vertices[2], fromIndex: vertices[1], target: &ntY);
        vectorBetweenVertices(toIndex: vertices[3], fromIndex: vertices[2], target: &ntZ);

        var x = ntX;
        var y = ntY;
        var z = ntZ;

        // This was generated using Mathematica
        let nx = x[3] * (y[2] * z[1] - y[1] * z[2])
            + x[2] * (y[1] * z[3] - y[3] * z[1])
            + x[1] * (y[3] * z[2] - y[2] * z[3]);
        let ny = x[3] * (y[0] * z[2] - y[2] * z[0])
            + x[2] * (y[3] * z[0] - y[0] * z[3])
            + x[0] * (y[2] * z[3] - y[3] * z[2]);
        let nz = x[3] * (y[1] * z[0] - y[0] * z[1])
            + x[1] * (y[0] * z[3] - y[3] * z[0])
            + x[0] * (y[3] * z[1] - y[1] * z[3]);
        let nw = x[2] * (y[0] * z[1] - y[1] * z[0])
            + x[1] * (y[2] * z[0] - y[0] * z[2])
            + x[0] * (y[1] * z[2] - y[2] * z[1]);

        let norm = sqrt(nx * nx + ny * ny + nz * nz + nw * nw);

        let f = 1.0 / norm;
        normal[0] = f * nx;
        normal[1] = f * ny;
        normal[2] = f * nz;
        normal[3] = f * nw;
    }

    /// Finds the normal vector nd.
    private func findNormalVectorND(vertices: [Int], normal: inout [Double]) {
        /* We need to solve the matrix A n = B where
         *  - A contains coordinates of vertices as columns
         *  - B is vector with all 1's. Really, it should be the distance of
         *      the plane from the origin, but - since we're not worried about that
         *      here and we will normalize the normal anyway - all 1's suffices.
         */
        var norm = 0.0;

        // Solve determinants by replacing x-th column by all 1.
        for x in 0..<dimension {
            for i in 0..<dimension {
                let offset = vertices[i] * dimension;
                for j in 0..<dimension {
                    // maybe I got the i/j mixed up here regarding the representation Math.net uses...
                    // ...but it does not matter since Det(A) = Det(Transpose(A)).
                    nDMatrix[dimension * i + j] = j == x ? 1.0 : positionData[offset + j];
                }
            }
            MathHelper.LUFactor(data: &nDMatrix, order: dimension, ipiv: &matrixPivots, vecLUcolj: &nDNormalHelperVector);
            var coord = 1.0;
            for i in 0..<dimension {
                if matrixPivots[i] != i {
                    coord *= -nDMatrix[dimension * i + i]; // the determinant sign changes on row swap.
                } else {
                    coord *= nDMatrix[dimension * i + i];
                }
            }
            normal[x] = coord;
            norm += coord * coord;
        }

        // Normalize the result
        let f = 1.0 / sqrt(norm);
        for i in 0..<normal.count {
            normal[i] *= f
        }
    }

    /// Gets the simplex volume. Prior to having enough edge vectors, the method pads the remaining with all
    /// "other numbers". So, yes, this method is not really finding the volume. But a relative volume-like measure. It
    /// uses the magnitude of the determinant as the volume stand-in following the Cayley-Menger theorem.
    internal func getSimplexVolume(edgeVectors: [[Double]], lastIndex: Int, bigNumber: Double) -> Double {
        var A = [Double](repeating: 0, count: dimension * dimension)
        var index = 0;
        for i in 0..<dimension {
            for j in 0..<dimension {
                if i <= lastIndex {
                    A[index] = edgeVectors[i][j];
                    index += 1
                } else {
                    A[index] = ( Double(truncating: pow(-1, index) as NSNumber) * Double(index)) / bigNumber
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
        return abs(determinantDestructive(&A));
    }

    /// Determinants the destructive.
    static var iPiv = [Int](repeating: 0, count: 5)
    static var helper = [Double](repeating: 0, count:5)

    private func determinantDestructive(_ A: inout [Double]) -> Double{
        switch dimension {
        case 0:
            return 0
        case 1:
            return A[0]
        case 2:
            return A[0] * A[3] - A[1] * A[2];
        case 3:
            return A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7]
                - A[0] * A[5] * A[7] - A[1] * A[3] * A[8] - A[2] * A[4] * A[6]
        default:

            MathHelper.LUFactor(data: &A, order: dimension, ipiv: &MathHelper.iPiv, vecLUcolj: &MathHelper.helper);
            var det = 1.0;
            for i in 0..<dimension {
                det *= A[dimension * i + i];
                if (MathHelper.iPiv[i] != i) {
                    det *= -1; // the determinant sign changes on row swap.
                }
            }
            return det;

        }
    }

    private static func LUFactor(data: inout [Double], order: Int, ipiv: inout [Int], vecLUcolj: inout [Double]) {
        // Initialize the pivot matrix to the identity permutation.
        for i in 0..<order {
            ipiv[i] = i
        }

        // Outer loop.
        for j in 0..<order {
            let indexj = j * order;
            let indexjj = indexj + j;

            // Make a copy of the j-th column to localize references.
            for i in 0..<order {
                vecLUcolj[i] = data[indexj + i];
            }

            // Apply previous transformations.
            for i in 0..<order {
                // Most of the time is spent in the following dot product.
                let kmax = min(i, j)
                var s = 0.0;
                for k in 0..<kmax {
                    s += data[k * order + i] * vecLUcolj[k];
                }
                vecLUcolj[i] -= s
                data[indexj + i] = vecLUcolj[i]
            }

            // Find pivot and exchange if necessary.
            var p = j;
            for i in (j + 1)..<order {
                if abs(vecLUcolj[i]) > abs(vecLUcolj[p]) {
                    p = i;
                }
            }

            if p != j {
                for k in 0..<order {
                    let indexk = k * order;
                    let indexkp = indexk + p;
                    let indexkj = indexk + j;
                    let temp = data[indexkp];
                    data[indexkp] = data[indexkj];
                    data[indexkj] = temp;
                }

                ipiv[j] = p;
            }

            // Compute multipliers.
            if j < order && data[indexjj] != 0.0 {
                for i in (j + 1)..<order {
                    data[indexj + i] /= data[indexjj];
                }
            }
        }
    }
}
