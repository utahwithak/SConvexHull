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
public class ConvexHull<TVertex: Vertex, TFace: ConvexFace> {

    /// Can only be created using a factory method.
    internal init() {
        points = []
        faces = []
    }
    internal init(points: [Vertex], faces: [TFace]) {
        self.points = points
        self.faces = faces
    }

    /// Points of the convex hull.
    public let points: [Vertex]

    /// Faces of the convex hull.
    public let faces: [TFace]

    /// Creates the convex hull
    public static func create<V: Vertex,F: ConvexFace>(with data: [V], planeDistanceTolerance tolerance: Double = Constants.defaultPlaneDistanceTolerance) -> ConvexHull<V,F> {
        return ConvexHullAlgorithm.getConvexHull(with: data, planeDistanceTolerance: tolerance)
    }

    public static func create<V: Vertex>(from data: [V], planeDistanceTolerance tolerance: Double = Constants.defaultPlaneDistanceTolerance) -> ConvexHull<V, DefaultConvexFace<V>> {
        return ConvexHull<V, DefaultConvexFace<TVertex>>.create(with: data, planeDistanceTolerance: tolerance)
    }


    /// Creates a convex hull of the input data.
    public static func create(from data: [[Double]], planeDistanceTolerance tolerance: Double = Constants.defaultPlaneDistanceTolerance) -> ConvexHull<DefaultVertex, DefaultConvexFace<DefaultVertex>> {
        let points = data.map{ DefaultVertex(position: $0)}
        return ConvexHull<DefaultVertex, DefaultConvexFace<DefaultVertex>>.create(with: points, planeDistanceTolerance: tolerance)
    }



}
