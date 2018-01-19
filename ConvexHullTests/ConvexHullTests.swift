//
//  ConvexHullTests.swift
//  ConvexHullTests
//
//  Created by Carl Wieland on 1/18/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import XCTest
@testable import ConvexHull

class ConvexHullTests: XCTestCase {
    
    override func setUp() {
        super.setUp()
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }
    
    override func tearDown() {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
        super.tearDown()
    }
    
    func testExample() {

        let verts = [[0.0,0.0,0.0],
                      [1.0,1.0,1.0],
                      [0.0,1.0,1.0],
                      [0.0,0.0,1.0],
                      [1.0,0.0,1.0],
                      [0.0,1.0,0.0],
                      [1.0,0.0,0.0],
                      [1.0,1.0,0.0]]
        let hull = ConvexHull<DefaultVertex, DefaultConvexFace<DefaultVertex>>.create(from: verts)
        XCTAssert(hull.faces.count == 12)
        // Use XCTAssert and related functions to verify your tests produce the correct results.
    }
    
    func testPerformanceExample() {
        // This is an example of a performance test case.
        self.measure {
            // Put the code you want to measure the time of here.
        }
    }
    
}
