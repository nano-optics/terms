library(ggplot2)
library(dplyr)
library(glue)
library(purrr)
library(tidyr)



ge <- get_geometry('input')
ge[,2:4] <-  as.matrix(ge[,2:4]) %*% terms::rotation_euler_active(-pi/4,0,0)

tpl_head <- '<?xml version="1.0" encoding="utf-8"?>
<scene version="0.5.0">
	<integrator  type="irrcache">
	  <integrator  type="photonmapper"/>
	</integrator>
'

tpl <- '
	<shape type="sphere">
	<point name="center" x="{x}" y="{y}" z="{z}"/> <float name="radius" value="{r}"/>

	<bsdf  type="roughconductor">
		<string  name="material"  value="{tag}"/>
		<float  name="alpha"  value="0.2"/>
	</bsdf>

	</shape>
'

tpl_foot <- '
	<sensor type="perspective" id="Camera-camera">
		<string name="fovAxis" value="smaller"/>
		<float name="focusDistance" value="50.0"/>
		<float name="fov" value="5"/>
		<transform name="toWorld">
			<lookAt target="0, 0, 0" origin="5000, 0, 0" up="0, 0, 10"/>
		</transform>

		<sampler type="ldsampler">
			<integer name="sampleCount" value="64"/>
		</sampler>

		<film type="hdrfilm" id="film">
			<integer name="width" value="400"/>
			<integer name="height" value="600"/>
			<string name="pixelFormat" value="rgb"/>
			<boolean name="banner" value="false"/>

			<rfilter type="gaussian"/>
		</film>
	</sensor>

	<emitter type="envmap" id="Area_002-light">
		<string name="filename" value="cloudy.hdr"/>
		<transform name="toWorld">
			<rotate y="1" angle="-180"/>
			<matrix value="-0.224951 -0.000001 -0.974370 0.000000 -0.974370 0.000000 0.224951 0.000000 0.000000 1.000000 -0.000001 8.870000 0.000000 0.000000 0.000000 1.000000 "/>
		</transform>
		<float name="scale" value="5"/>
	</emitter>

	<bsdf type="diffuse" id="__diffmat">
		<rgb name="reflectance" value="0.3 0.3 0.3"/>
	</bsdf>

	<texture type="checkerboard" id="__planetex">
		<rgb name="color0" value="#e0ebff"/>
		<rgb name="color1" value="#e0ebff"/>
		<float name="uscale" value="100.0"/>
		<float name="vscale" value="100.0"/>
		<float name="uoffset" value="0.0"/>
		<float name="voffset" value="100.0"/>
	</texture>

	<bsdf type="diffuse" id="__planemat">
		<ref name="reflectance" id="__planetex"/>
	</bsdf>


	<shape type="serialized" id="Plane-mesh_0">
		<string name="filename" value="matpreview.serialized"/>
		<integer name="shapeIndex" value="0"/>
		<transform name="toWorld">
			<rotate z="1" angle="90"/>
			<matrix value="1000 -1000 0 -1000 1000 1000 0 1000 0 0 1000 -100 0 0 0 1"/>
				<translate z="-500"/>
		</transform>

		<ref name="bsdf" id="__planemat"/>
	</shape>
</scene>
'

out <- 'vis_helix.xml'
cat(tpl_head, '\n', file = out)
cat(glue_data(tpl, .x = ge), '\n', file = out, append = T)
cat(tpl_foot, '\n', file = out, append = T)

