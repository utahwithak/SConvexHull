//
//  Vector3.swift
//  ConvexHull
//
//  Created by Carl Wieland on 1/19/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation

public struct Vector3 {
    public var position: [Double]
    public var x: Double {
        get {
            return position[0]
        }
        set {
            position[0] = newValue
        }
    }
    public var y: Double {
        get {
            return position[1]
        }
        set {
            position[1] = newValue
        }
    }
    public var z: Double {
        get {
            return position[2]
        }
        set {
            position[2] = newValue
        }
    }
    public init(x: Int, y:Int, z: Int) {
        position = [Double(x),Double(y),Double(z)]
    }

    public init(x: Double, y: Double, z: Double) {
        position = [x,y,z]
    }
    
    public init(x: CGFloat, y: CGFloat, z: CGFloat) {
        position = [Double(x),Double(y),Double(z)]
    }
    public init(positions: [Double]) {
        self.position = positions
    }

    subscript(index: Int) -> Double {
        get {
            return position[index]
        }
        set {
            position[index] = newValue
        }
    }
}
