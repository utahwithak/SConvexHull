//
//  Vector3.swift
//  ConvexHull
//
//  Created by Carl Wieland on 1/19/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation

public struct Vector3 {
    public static let zero = Vector3(x: 0, y: 0, z: 0)
    public var x: Double
    public var y: Double
    public var z: Double

    public init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
    
    public init(x: CGFloat, y: CGFloat, z: CGFloat) {
        self.init(x: Double(x), y: Double(y),z: Double(z))
    }

    public init(x: Int, y:Int, z: Int) {
        self.init(x: Double(x), y: Double(y),z: Double(z))
    }

}

func *(lhs: Double, rhs: Vector3) -> Vector3 {
    return Vector3(x: lhs * rhs.x, y: lhs * rhs.y, z: lhs * rhs.z)
}

func *(lhs: Vector3, rhs: Double) -> Vector3 {
    return Vector3(x: lhs.x * rhs, y: lhs.y * rhs, z: lhs.z * rhs)
}

func *=(lhs: inout Vector3, rhs: Double) {
    lhs.x *= rhs
    lhs.y *= rhs
    lhs.z *= rhs
}

func +=(lhs: inout Vector3, rhs: Vector3) {
    lhs.x += rhs.x
    lhs.y += rhs.y
    lhs.z += rhs.z
}

func -(lhs: Vector3, rhs: Vector3) -> Vector3{
    return Vector3(x: lhs.x - rhs.x, y: lhs.y - rhs.y, z: lhs.z - rhs.z)
}

func +(lhs: Vector3, rhs: Vector3) -> Vector3{
    return Vector3(x: lhs.x + rhs.x, y: lhs.y + rhs.y, z: lhs.z + rhs.z)
}
