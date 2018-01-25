//
//  ConvexHull.swift
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

/// Representation of a convex hull.
public struct ConvexHull {
    /// The default plane distance tolerance
    public static let defaultPlaneDistanceTolerance: Double = 1e-10
    internal init(points: [Vector3], faces: [ConvexFace]) {
        self.points = points
        self.faces = faces
    }

    /// Points of the convex hull.
    public let points: [Vector3]

    /// Faces of the convex hull.
    public let faces: [ConvexFace]

    /// Creates the convex hull
    public static func create(with data: [Vector3], planeDistanceTolerance tolerance: Double = ConvexHull.defaultPlaneDistanceTolerance) -> ConvexHull{
        return ConvexHullAlgorithm.getConvexHull(with: data, planeDistanceTolerance: tolerance)
    }

}
