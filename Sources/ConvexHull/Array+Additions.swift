//
//  Array+Additions.swift
//  ConvexHull
//
//  Created by Carl Wieland on 1/19/18.


import Foundation

extension Array {

    /// Removes all elements where `predicate` returns true
    ///
    /// - Parameter predicate: Closure to test elements against. Removes item if it `true`.
    mutating func removeWhere(_ predicate: (Element) -> Bool) {
        for (index, element) in self.enumerated().reversed() {
            if predicate(element) {
                remove(at: index)
            }
        }
    }
}
