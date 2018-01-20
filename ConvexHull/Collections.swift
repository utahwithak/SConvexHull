//
//  Collections.swift
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


/// A more lightweight alternative to List of T.
/// On clear, only resets the count and does not clear the references
/// =&gt; this works because of the ObjectManager.
/// Includes a stack functionality.

internal class SimpleList<T> {

    /// The count
    public var count: Int {
        return items.count
    }

    /// The items
    private var items = [T]()

    /// Get the i-th element.
    subscript (index: Int) -> T {
        get {
            return items[index]
        }
        set {
            items[index] = newValue
        }
    }

    /// Adds a vertex to the buffer.
    public func append(_ item: T){
        items.append(item)
    }

    /// Pops the last value from the list.
    public func pop() -> T? {
        return items.popLast()
    }

    /// Sets the Count to 0, otherwise does nothing.
    public func clear() {
        items.removeAll(keepingCapacity: true)
    }
}

/// Class IndexBuffer.
/// A fancy name for a list of integers.
internal class IndexBuffer : SimpleList<Int> { }

/// A priority based linked list.
internal final class FaceList {

    /// The last
    private var last: ConvexFaceInternal?;


    /// Get the first element.

    public var first: ConvexFaceInternal?


    /// Adds the element to the beginning.

    private func addFirst(_ face: ConvexFaceInternal) {
        face.isInList = true;
        first?.previous = face;
        face.next = first
        first = face;
    }

    /// Adds a face to the list.
    public func append(_ face: ConvexFaceInternal) {
        if face.isInList {

            if first!.verticesBeyond.count < face.verticesBeyond.count {
                remove(face);
                addFirst(face);
            }
            return;
        }


        face.isInList = true;

        if let first = first, first.verticesBeyond.count < face.verticesBeyond.count {
            first.previous = face;
            face.next = first;
            self.first = face;
        } else {
            if let last = last {
                last.next = face
            }
            face.previous = last
            last = face
            if first == nil {
                first = face;
            }
        }
    }


    /// Removes the element from the list.
    public func remove(_ face: ConvexFaceInternal) {

        guard face.isInList else {
            return
        }

        face.isInList = false;

        if let prev = face.previous {
            prev.next = face.next;
        } else if face.previous == nil {
            first = face.next;
        }

        if let next = face.next {
            next.previous = face.previous
        } else if face.next == nil {
            last = face.previous
        }

        face.next = nil
        face.previous = nil
    }
}

internal final class ConnectorList {

    private var last: FaceConnector?

    public var first: FaceConnector?

    public func append(_ element: FaceConnector) {
        if let last = last {
            last.next = element;
        }

        element.previous = last;
        last = element;
        if first == nil {
            first = element
        }
    }

    public func remove( connector: FaceConnector) {

        if let prev = connector.previous {
            prev.next = connector.next
        } else if connector.previous == nil {
            first = connector.next
        }

        if let next = connector.next {
            next.previous = connector.previous;
        } else if connector.next == nil {
            last = connector.previous
        }

        connector.next = nil
        connector.previous = nil
    }
}
