"""European city lookup using Open-Meteo / GeoNames, with a local cache."""

import hashlib
import json
import math
import re
import unicodedata
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError


GEOCODING_URL = "https://geocoding-api.open-meteo.com/v1/search"
EUROPEAN_ISLAND_TIMEZONES = {
    "Atlantic/Reykjavik", "Atlantic/Faroe", "Atlantic/Faeroe",
    "Atlantic/Azores", "Atlantic/Madeira", "Atlantic/Canary",
    "Arctic/Longyearbyen", "Asia/Nicosia", "Asia/Famagusta",
}


def city_slug(city):
    """Keep Unicode city names readable and remove unsafe path characters."""
    name = unicodedata.normalize("NFKC", city).casefold().strip()
    name = re.sub(r"[^\w-]+", "_", name, flags=re.UNICODE).strip("_-")
    if not name:
        raise ValueError("City must contain letters or numbers")
    # Prefix reserved Windows device names, including when a country is absent.
    if re.fullmatch(r"con|prn|aux|nul|com[1-9]|lpt[1-9]", name):
        name = "city_" + name
    return name[:100]


def validate_location(latitude, longitude, timezone):
    if not (math.isfinite(latitude) and math.isfinite(longitude)
            and -90 <= latitude <= 90 and -180 <= longitude <= 180):
        raise ValueError("Invalid latitude/longitude")
    if not (timezone.startswith("Europe/") or timezone in EUROPEAN_ISLAND_TIMEZONES):
        raise ValueError("Choose a European city and its local IANA timezone")
    try:
        ZoneInfo(timezone)
    except ZoneInfoNotFoundError as exc:
        raise ValueError(f"Unknown timezone {timezone!r}; install the tzdata package") from exc


def lookup_city(city, country_code=None, cache_dir="pvgis_cache", refresh=False):
    """Resolve a populated place; qualifiers and ISO country codes narrow searches.

    If names are ambiguous, choose the most populous exact match in Europe.
    The resolved place is always reported by the workflow and stored in metadata.
    """
    city = city.strip()
    if len(city.split(",", 1)[0].strip()) < 2:
        raise ValueError("Provide a city name of at least two characters")
    params = {"name": city, "count": 100, "language": "en", "format": "json"}
    if country_code:
        country_code = country_code.strip().upper()
        if not re.fullmatch(r"[A-Z]{2}", country_code):
            raise ValueError("country_code must be a two-letter ISO code, e.g. NO")
        params["countryCode"] = country_code
    key = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:16]
    cache = Path(cache_dir) / "locations" / f"{city_slug(city)}_{key}.json"
    if cache.exists() and not refresh:
        payload = json.loads(cache.read_text(encoding="utf-8"))
    else:
        request = Request(GEOCODING_URL + "?" + urlencode(params),
                          headers={"User-Agent": "pv-battery-sizing/1.0"})
        try:
            with urlopen(request, timeout=30) as response:
                payload = json.load(response)
        except (URLError, TimeoutError, OSError) as exc:
            raise RuntimeError(
                f"Could not resolve city {city!r}. Check connectivity, use a cached "
                "location, or provide latitude, longitude and timezone explicitly."
            ) from exc
        if payload.get("error"):
            raise ValueError(f"City lookup failed: {payload.get('reason')}")
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    matches = []
    for place in payload.get("results", []):
        if not place.get("feature_code", "").startswith("P"):
            continue
        if country_code and place.get("country_code") != country_code:
            continue
        try:
            validate_location(place["latitude"], place["longitude"], place["timezone"])
        except (KeyError, ValueError):
            continue
        matches.append(place)
    if not matches:
        raise ValueError(f"No European city found for {city!r}. Try a local spelling, "
                         "a country code, or 'City, Region'.")

    def normalized(value):
        return "".join(ch for ch in unicodedata.normalize("NFKD", value.casefold())
                       if not unicodedata.combining(ch))

    query = normalized(city.split(",", 1)[0].strip())
    exact = [p for p in matches if normalized(p["name"]) == query]
    place = max(exact or matches, key=lambda p: p.get("population", 0))
    return {
        "latitude_deg": float(place["latitude"]),
        "longitude_deg": float(place["longitude"]),
        "local_timezone": place["timezone"],
        "country_code": place.get("country_code", ""),
        "resolved_city": place["name"],
        "location_description": ", ".join(str(place[k]) for k in ("name", "admin1", "country")
                                           if place.get(k)),
        "location_source": f"Open-Meteo / GeoNames (ID {place['id']})",
    }
