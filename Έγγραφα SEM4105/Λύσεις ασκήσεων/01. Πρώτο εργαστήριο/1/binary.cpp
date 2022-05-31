#include <iostream>
#include <array>
#include <algorithm>
#include <iterator>

int main()
{
    std::cout << u8"Δώσε μη αρνητικό ακέραιο: ";

    int k;
    std::cin >> k;

    std::array<int,32> b;

    std::generate(b.rbegin(), b.rend(),
		  [&k] ()
		      {
			  auto x = k%2;
			  k/=2;
			  return x;
		      }
	);

    auto beg = std::find(b.cbegin(), b.cend(), 1);

    std::ostream_iterator<int> out{std::cout, ""};
    std::copy(beg, b.cend(), out);
    std::cout << '\n';
}
