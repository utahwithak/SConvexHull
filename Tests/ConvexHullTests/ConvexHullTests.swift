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
        
        let verts = [Vector3(0.0,0.0,0.0),
                     Vector3(1.0,1.0,1.0),
                     Vector3(0.0,1.0,1.0),
                     Vector3(0.0,0.0,1.0),
                     Vector3(1.0,0.0,1.0),
                     Vector3(0.0,1.0,0.0),
                     Vector3(1.0,0.0,0.0),
                     Vector3(1.0,1.0,0.0)]
        var hull: ConvexHull? = ConvexHull.create(with: verts)
        
        XCTAssert(hull!.faces.count == 12)
        
        hull = nil
        
        // Expected that memory doesn't leak
        XCTAssertNotNil(hull)
    }
    
    func testPerformanceExample() {
        // This is an example of a performance test case.
        self.measure {
            // Put the code you want to measure the time of here.
        }
    }
    
}
