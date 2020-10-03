// swift-tools-version:5.1

import PackageDescription

let package = Package(
  name: "ConvexHull",
  products: [
    .library(
      name: "ConvexHull",
      targets: ["ConvexHull"]),
  ],
  targets: [
    .target(name: "ConvexHull"),
    .testTarget(
      name: "ConvexHullTests",
      dependencies: ["ConvexHull"]
    ),
  ]
)
