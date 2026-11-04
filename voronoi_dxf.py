
import csv
import ezdxf
import numpy as np

from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
#datos de entrada
DXF_PUNTOS = "puntos.dxf"
DXF_REGION = "region_exterior.dxf"
DXF_ZONA1 = "zona_influencia_1.dxf"
DXF_ZONA2 = "zona_influencia_2.dxf"
#salida
DXF_SALIDA = "resultado_voronoi.dxf"
CSV_SALIDA = "areas_por_punto.csv"
#leer los poligonos
def leer_poligonos(ruta):
    doc = ezdxf.readfile(ruta)
    msp = doc.modelspace()
    poligonos = []

    for entidad in msp:
        tipo = entidad.dxftype()

        if tipo == "LWPOLYLINE":
            if not entidad.closed:
                continue

            coords = [
                (p[0], p[1])
                for p in entidad.get_points()
            ]

        elif tipo == "POLYLINE":
            if not entidad.is_closed:
                continue

            coords = [
                (v.dxf.location.x, v.dxf.location.y)
                for v in entidad.vertices
            ]

        else:
            continue

        if len(coords) < 3:
            continue

        poligono = Polygon(coords)

        if not poligono.is_valid:
            poligono = poligono.buffer(0)

        if not poligono.is_empty and poligono.area > 0:
            poligonos.append(poligono)

    if not poligonos:
        raise ValueError(
            f"No se encontraron polilineas cerradas en {ruta}"
        )

    return unary_union(poligonos)


# LECTURA DE PUNTOS DXF

def leer_puntos(ruta):
    doc = ezdxf.readfile(ruta)
    puntos = []

    for entidad in doc.modelspace().query("POINT"):
        x = entidad.dxf.location.x
        y = entidad.dxf.location.y

        puntos.append((x, y))

    if len(puntos) < 2:
        raise ValueError(
            "Se necesitan al menos dos puntos distintos."
        )

    unicos = list(dict.fromkeys(puntos))

    if len(unicos) < 2:
        raise ValueError(
            "No hay suficientes puntos distintos."
        )

    return unicos

# Voronoi regiones delimitadas 

def voronoi_finito(puntos, radio):
    vor = Voronoi(np.asarray(puntos))
    regiones = []
    vertices = vor.vertices.tolist()

    centro = np.mean(vor.points, axis=0)
    crestas = {}

    for (i, j), (v1, v2) in zip(
        vor.ridge_points, vor.ridge_vertices
    ):
        crestas.setdefault(i, []).append((j, v1, v2))
        crestas.setdefault(j, []).append((i, v1, v2))

    for i, region_id in enumerate(vor.point_region):
        indices = vor.regions[region_id]

        if indices and all(v >= 0 for v in indices):
            regiones.append(
                Polygon([vertices[v] for v in indices])
            )
            continue

        nuevos_indices = [
            v for v in indices if v >= 0
        ] if indices else []

        for j, v1, v2 in crestas.get(i, []):
            if v1 >= 0 and v2 >= 0:
                continue

            v = v1 if v1 >= 0 else v2

            direccion = vor.points[j] - vor.points[i]
            norma = np.linalg.norm(direccion)

            if norma == 0:
                continue

            direccion = direccion / norma
            normal = np.array(
                [-direccion[1], direccion[0]]
            )

            punto_medio = (
                vor.points[i] + vor.points[j]
            ) / 2

            if np.dot(
                punto_medio - centro, normal
            ) < 0:
                normal = -normal

            extremo = vor.vertices[v] + normal * radio
            vertices.append(extremo.tolist())
            nuevos_indices.append(len(vertices) - 1)

        coords = [vertices[v] for v in nuevos_indices]

        if len(coords) < 3:
            raise ValueError(
                "No se pudo cerrar una celda de Voronoi."
            )

        centroide = np.mean(coords, axis=0)

        coords.sort(
            key=lambda p: np.arctan2(
                p[1] - centroide[1],
                p[0] - centroide[0]
            )
        )

        regiones.append(Polygon(coords))

    return regiones


# Exportar la geometria a formato dxf

def escribir_geometria(msp, geometria, capa):
    if geometria.is_empty:
        return

    if isinstance(geometria, Polygon):
        poligonos = [geometria]
    elif isinstance(geometria, MultiPolygon):
        poligonos = list(geometria.geoms)
    else:
        poligonos = [
            g for g in getattr(geometria, "geoms", [])
            if isinstance(g, Polygon)
        ]

    for poligono in poligonos:
        if poligono.area <= 0:
            continue

        msp.add_lwpolyline(
            list(poligono.exterior.coords)[:-1],
            close=True,
            dxfattribs={"layer": capa}
        )

#funcion centroide de las figuras

def escribir_centroide(msp, geometria, capa, identificador):
    if geometria.is_empty or geometria.area <= 0:
        return None

    centroide = geometria.centroid

    x = centroide.x
    y = centroide.y

    # Dibujar el centroide como un punto
    msp.add_point(
        (x, y),
        dxfattribs={"layer": capa}
    )

    # Agregar una etiqueta con el identificador
    msp.add_text(
        f"{identificador}: ({x:.3f}, {y:.3f})",
        dxfattribs={
            "layer": capa,
            "height": 0.5
        }
    )

    return x, y

#funcion main

def main():
    region = leer_poligonos(DXF_REGION)
    zona1 = leer_poligonos(DXF_ZONA1)
    zona2 = leer_poligonos(DXF_ZONA2)
    puntos = leer_puntos(DXF_PUNTOS)

    # limitar las zonas
    zona1 = zona1.intersection(region)
    zona2 = zona2.intersection(region)

    
    minx, miny, maxx, maxy = region.bounds
    radio = 100 * max(maxx - minx, maxy - miny, 1)

    celdas = voronoi_finito(puntos, radio)

    salida = ezdxf.new("R2010")
    msp = salida.modelspace()

    for capa in [
        "REGION",
        "VORONOI",
        "INTERSECCION_Z1",
        "INTERSECCION_Z2",
        "CENTROIDES_VORONOI",
        "CENTROIDES_Z1",
        "CENTROIDES_Z2"
    ]:
        salida.layers.new(capa)

    escribir_geometria(msp, region, "REGION")

    resultados = []

    for i, (xy, celda) in enumerate(zip(puntos, celdas), start=1):
        # Recorte de la celda por la region exterior.
        celda = celda.intersection(region)

        # Interseccion de cada celda con cada zona.
        parte1 = celda.intersection(zona1)
        parte2 = celda.intersection(zona2)

        area_total = celda.area
        area1 = parte1.area
        area2 = parte2.area

        

        escribir_geometria(msp, celda, "VORONOI")
        escribir_geometria(
            msp, parte1, "INTERSECCION_Z1"
        )
        escribir_geometria(
            msp, parte2, "INTERSECCION_Z2"
        )
        centro_v = escribir_centroide(
            msp, celda, "CENTROIDES_VORONOI", f"V{i}"
        )

        centro_z1 = escribir_centroide(
            msp, parte1, "CENTROIDES_Z1", f"Z1_{i}"
        )

        centro_z2 = escribir_centroide(
            msp, parte2, "CENTROIDES_Z2", f"Z2_{i}"
        )

        resultados.append({
            "id_punto": i,
            "x": xy[0],
            "y": xy[1],
            "area_voronoi": area_total,
            "area_zona_1": area1,
            "area_zona_2": area2,
            "centroide_voronoi_x": centro_v[0] if centro_v else None,
            "centroide_voronoi_y": centro_v[1] if centro_v else None,

            "centroide_zona1_x": centro_z1[0] if centro_z1 else None,
            "centroide_zona1_y": centro_z1[1] if centro_z1 else None,

            "centroide_zona2_x": centro_z2[0] if centro_z2 else None,
            "centroide_zona2_y": centro_z2[1] if centro_z2 else None,
        })

    salida.saveas(DXF_SALIDA)

    with open(
        CSV_SALIDA, "w", newline="", encoding="utf-8-sig"
    ) as archivo:
        campos = [
            "id_punto", "x", "y",
            "area_voronoi", "area_zona_1", "area_zona_2",
            "centroide_voronoi_x", "centroide_voronoi_y",
            "centroide_zona1_x", "centroide_zona1_y",
            "centroide_zona2_x", "centroide_zona2_y"
        ]

        escritor = csv.DictWriter(
            archivo, fieldnames=campos
        )
        escritor.writeheader()
        escritor.writerows(resultados)
    
        print("Proceso terminado.")
        print(f"DXF: {DXF_SALIDA}")
        print(f"CSV: {CSV_SALIDA}")

        for r in resultados:
            print(
                f"Punto {r['id_punto']}: "
                f"Voronoi={r['area_voronoi']:.3f}, "
                f"Zona 1={r['area_zona_1']:.3f}, "
                f"Zona 2={r['area_zona_2']:.3f}"
            )


if __name__ == "__main__":
    main()


