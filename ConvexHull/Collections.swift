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
