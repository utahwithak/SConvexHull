//
//  DataModel.swift
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

internal struct DeferredFace {

    /// The faces.
    let face: ConvexFaceInternal
    let pivot: ConvexFaceInternal
    let oldFace: ConvexFaceInternal


    /// The indices.
    let faceIndex: Int
    let pivotIndex: Int
}


/// A helper class used to connect faces.
internal struct FaceConnector {

    /// The edge to be connected.
    let edgeIndex: Int

    /// The face.
    let face: ConvexFaceInternal

    /// The hash code computed from indices.
    let hashCode: UInt64

    /// The vertex indices.
    let v0: Int
    let v1: Int

    init(face: ConvexFaceInternal, edgeIndex: Int) {
        self.face = face;
        self.edgeIndex = edgeIndex;

        switch edgeIndex {
        case 0:
            v0 = face.vertices[1]
            v1 = face.vertices[2]
        case 1:
            v0 = face.vertices[0]
            v1 = face.vertices[2]
        default:
            v0 = face.vertices[0]
            v1 = face.vertices[1]
        }

        var hashCode = UInt64(23)
        hashCode = hashCode &+ (31 &* hashCode &+ UInt64(v0))
        hashCode = hashCode &+ (31 &* hashCode &+ UInt64(v1))

        self.hashCode = hashCode
    }

}


/// This internal class manages the faces of the convex hull. It is a
/// separate class from the desired user class.
internal final class ConvexFaceInternal
{
    /// Gets or sets the adjacent face data.
    public var adjacentFaces = [0,0,0]

    /// The furthest vertex.
    public var furthestVertex = 0

    /// Index of the face inside the pool.
    public var index = 0

    /// Is it present in the list.
    public var isInList = false

    /// Is the normal flipped?
    public var isNormalFlipped = false

//    /// Next node in the list.
//    public var next: ConvexFaceInternal?

    /// Gets or sets the normal vector.
    public var normal = Vector3.zero

    /// Face plane constant element.
    public var offset: Double = 0

//    /// Prev node in the list.
//    public weak var previous: ConvexFaceInternal?

    /// Used to traverse affected faces and create the Delaunay representation.
    public var tag = 0

    /// Gets or sets the vertices.
    public var vertices = [0,0,0]

    /// Gets or sets the vertices beyond.
    public var verticesBeyond = [Int]()

    public init(index: Int) {
        self.index = index
    }

    func reset() {
        for i in 0..<adjacentFaces.count {
            adjacentFaces[i] = -1
        }
    }
}
