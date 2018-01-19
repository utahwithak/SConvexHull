//
//  ObjectManager.swift
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

/// A helper class for object allocation/storage.
/// This helps the GC a lot as it prevents the creation of about 75% of
/// new face objects (in the case of ConvexFaceInternal). In the case of
/// FaceConnectors and DefferedFaces, the difference is even higher (in most
/// cases O(1) vs O(number of created faces)).
internal class ObjectManager {

    private let dimension: Int

    private var connectorStack: FaceConnector?

    private let deferredFaceStack: SimpleList<DeferredFace>

    private let emptyBufferStack: SimpleList<IndexBuffer>

    private let facePool: SimpleList<ConvexFaceInternal>

    private let freeFaceIndices: IndexBuffer

    init(dimension: Int, facePool: SimpleList<ConvexFaceInternal>) {
        self.dimension = dimension
        self.facePool = facePool;
        freeFaceIndices = IndexBuffer()
        emptyBufferStack = SimpleList<IndexBuffer>()
        deferredFaceStack = SimpleList<DeferredFace>()
    }

    /// Return the face to the pool for later use.
    internal func depositFace(at faceIndex: Int) {
        let face = facePool[faceIndex]
        face.reset()
        freeFaceIndices.push(faceIndex);
    }


    /// Create a new face and put it in the pool.
    private func createFace() -> Int {
        let index = facePool.count
        let face = ConvexFaceInternal(dimension: dimension, index: index, beyondList: getVertexBuffer());
        facePool.append(face)
        return index
    }

    public func getFace() -> Int {
        return freeFaceIndices.pop() ?? createFace()
    }

    /// Store a face connector in the "embedded" linked list.
    public func depositConnector(_ connector: FaceConnector) {
        if let stack = connectorStack {
            connector.next = stack
            connectorStack = connector;
        } else {
            connector.next = nil;
            connectorStack = connector;
        }


    }

    /// Get an unused face connector. If none is available, create it.
    public func getConnector() -> FaceConnector {

        if let ret = connectorStack {
            connectorStack = ret.next
            ret.next = nil
            return ret;
        } else {
            return FaceConnector(dimension: dimension)
        }
    }

    /// Deposit the index buffer.
    public func depositVertexBuffer(_ buffer: IndexBuffer) {
        buffer.clear();
        emptyBufferStack.append(buffer);
    }


    /// Get a store index buffer or create a new instance.
    public func getVertexBuffer() -> IndexBuffer {
        return emptyBufferStack.pop() ?? IndexBuffer()
    }

    /// Deposit the deferred face.
    public func depositDeferredFace(_ face: DeferredFace) {
        deferredFaceStack.append(face);
    }


    /// Get the deferred face.
    public func getDeferredFace() -> DeferredFace {
        return deferredFaceStack.pop() ?? DeferredFace()
    }

}
