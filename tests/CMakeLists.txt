CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ADD_EXECUTABLE (simpleTestTriangle simpleTestTriangle.cpp)
TARGET_LINK_LIBRARIES (simpleTestTriangle BasicCollision)

ADD_EXECUTABLE (simpleTestRectangle simpleTestRectangle.cpp)
TARGET_LINK_LIBRARIES (simpleTestRectangle BasicCollision)

ADD_EXECUTABLE (simpleTestBox simpleTestBox.cpp)
TARGET_LINK_LIBRARIES (simpleTestBox BasicCollision)

ADD_EXECUTABLE (simpleTestBoxes simpleTestBoxes.cpp)
TARGET_LINK_LIBRARIES (simpleTestBoxes BasicCollision)