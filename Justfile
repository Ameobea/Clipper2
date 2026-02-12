# Root-level tasks for this repo

# Build and run the C++ test suite (includes critical point tests)
cpp-tests:
  cmake --build CPP/build --target ClipperTests
  cd CPP/build && ./ClipperTests
