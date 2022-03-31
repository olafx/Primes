/*
A segmented sieve of Eratosthenes.

TODO
-   The starting from the square of the prime optimization is nice. Does a version of this also work for the other
    segments, i.e. the ones that don't start from the origin?
-   The size of prime0 can be halved by not considering even numbers. This is really quite strange. Would an
    optimization like this additionally work for odd numbers like 3? Need to better understand why and how this
    type of optimization works before implementing it.
*/

#include <iostream>
#include <cstdio>
#include <cmath>

int main(int argc, char **argv)
{
    /* primes up to and including n */ size_t n;
    if (argc != 2 || sscanf(argv[1], "%zu\n", &n) == 0 || n < 9)
        return EXIT_FAILURE;

    /* number of primes thusfar */ size_t N = 0;

    // The problem is split into segments [m - d, m) where d = floor(sqrt(n)). m is initially d. These segments have
    // length d, so m increases in steps of d. The last m must be at least n + 1 so that n is included. d will need to
    // be at least 3 for the rest to work, hence the requirement for n to be at least 9.
    /* segment size */ const size_t d = sqrt(n);

    // Find the primes in the first segment [0, d) using a sieve of Erathostenes.

    // Here prime0[i] is i being prime.
    /* segment sieve array */ auto *const prime0 = new bool[d];

    // Initial filling. Sieve starts at enumerations of the first prime: 2.
    // Since d is the largest prime we might possibly have to store in memory, we'd like to also include d in this first
    // segment (but not count it towards the prime count), so that we're done with all primes in memory right away.
    for (size_t i = 0; i < d - 1; i++)
        prime0[i] = true;

    // Sieve of Erathostenes.
    const size_t max1 = sqrt(d);
    for (size_t i = 2; i <= max1;)
    {   // Enumerate and mark the multiples of i. Can skip right away to i^2, the ones before are already done.
        for (size_t e = i * i; e <= d; e += i)
            prime0[e - 2] = false;
        // The next prime is the first unmarked one.
        for (i++; i <= d && !prime0[i - 2]; i++);
    }

    // Update number of primes thusfar. We don't want to include d in the count since it's not part of [0, d).
    for (size_t i = 2; i < d; i++)
        if (prime0[i - 2])
            N++;

    // The number of primes we'll be storing.
    const size_t N1 = N + (prime0[d - 2] ? 1 : 0);

    // Store these primes as actual numbers now so that prime1[i] is the i-th prime up to and including d.
    /* primes up to and including floor(sqrt(n)) */ auto *const prime1 = new size_t[N1];
    for (size_t i = 2, j = 0; i <= d; i++)
        if (prime0[i - 2])
            prime1[j++] = i;

    // Find the primes in the other segments [m - d, m) by enumerating the primes from the 1st interval through each.
    // This works because the first hole is definitely prime, and its square is not in the same interval [m - d, m), so
    // all holes will be prime. None of these have to be stored, we already have primes up to and including d stored.

    // Iterate over the other segments [m - d, m), up to potentially m = n; the last segment is an exception for after.
    size_t m;
    for (m = 2 * d; m <= n; m += d)
    {
        /* left extreme of segment */ const size_t l = m - d;
        const size_t max2 = sqrt(m);

        // Reset the sieve.
        for (size_t i = 0; i < d; i++)
            prime0[i] = true;

        // Sieve of Erathostenes.
        for (size_t i = 0; i < N1; i++)
        {   const size_t p = prime1[i];
            // Only considering primes up to and including floor(sqrt(m)).
            if (p > max2)
                break;
            // Calculate the smallest multiple of p in [m - d, m). This is ceil(l / p) * p.
            const size_t s = p * (l / p + (l % p != 0));
            // Enumerate and mark, starting at that smallest multiple.
            for (size_t e = s; e < m; e += p)
                prime0[e - l] = false;
        }

        // Update number of primes thusfar.
        for (size_t i = 0; i < d; i++)
            if (prime0[i])
                N++;
    }

    // Now do the last segment (the one that includes n) and stop at n. For this segment, m > n. Because m > n and
    // because it's the last segment, some optimizations can be made.
    {
        /* left extreme of segment */ const size_t l = m - d;
        // Don't need a max2.

        // Reset the sieve. Don't need to reset all of it.
        for (size_t i = 0; i <= n - l; i++)
            prime0[i] = true;

        // Sieve of Erathostenes.
        for (size_t i = 0; i < N1; i++)
        {   const size_t p = prime1[i];
            // Considering all primes up to and including d now, so no break.
            // Calculate the smallest multiple of p in [m - d, m). This is ceil(l / p) * p.
            const size_t s = p * (l / p + (l % p != 0));
            // Enumerate and mark up to and including n, starting at that smallest multiple.
            for (size_t e = s; e <= n; e += p)
                prime0[e - l] = false;
        }

        // Update number of primes thusfar. Stop at n.
        for (size_t i = 0; i <= n - l; i++)
            if (prime0[i])
                N++;
    }

    delete[] prime0;
    delete[] prime1;

    std::cout << N << '\n';
}
