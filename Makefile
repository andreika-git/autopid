pid_from_msl: pid_from_msl.cpp
	g++ -std=c++17 -fpermissive -O3 -static -o pid_from_msl.exe pid_from_msl.cpp googletest/gtest-all.cpp googletest/gmock-all.cpp -Igoogletest/googletest -Igoogletest/googletest/include -Igoogletest/googlemock/include
