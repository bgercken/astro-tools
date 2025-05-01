Computes dates/times when:

* The Moon is below the horizon
* Astronomical darkness is present (Sun below -18°)
* Milky Way core is above 15° elevation

*Setup:*

pip install skyfield numpy pandas

*Notes:*

1. It uses Skyfield and the JPL DE421 ephemeris.
2. For best results, run it regularly with updated dates.
3. It can be expanded to include weather when internet is available later.
