/* global define */
define('scripts/gl_helper', [], function() {
  'use strict';

  return class GlHelper {
    constructor(canvas) {
      this.canvas = canvas;

      this.gl = this.canvas.getContext('webgl2');
      this.gl.getExtension('EXT_color_buffer_float');
      this.gl.getExtension('OES_texture_float_linear');
      this.gl.getExtension('EXT_float_blend');
    }

    loadShader(type, source) {
      const gl = this.gl;

      const shader = gl.createShader(type);
      gl.shaderSource(shader, source);
      gl.compileShader(shader);

      if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert('An error occured compiling the shaders: ' + gl.getShaderInfoLog(shader));
        gl.deleteShader(shader);
        return null;
      }

      return shader;
    }

    loadFloatTexture(width, height, values) {
      const gl = this.gl;

      const texture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, texture);

      const level = 0;
      const internalFormat = gl.RGBA32F;
      const border = 0;
      const format = gl.RGBA;
      const type = gl.FLOAT;

      gl.texImage2D(gl.TEXTURE_2D, level, internalFormat, width, height, border, format, type, values);

      gl.texParameterf(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameterf(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameterf(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameterf(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

      gl.bindTexture(gl.TEXTURE_2D, null);

      return texture;
    }

    loadUintTexture(width, height, values) {
      const gl = this.gl;

      const texture = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, texture);

      const level = 0;
      const internalFormat = gl.RGBA32UI;
      const border = 0;
      const format = gl.RGBA_INTEGER;
      const type = gl.UNSIGNED_INT;

      gl.texImage2D(gl.TEXTURE_2D, level, internalFormat, width, height, border, format, type, values);

      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

      gl.bindTexture(gl.TEXTURE_2D, null);

      return texture;
    }

    loadShaderProgram(vertexShaderSource, fragmentShaderSource) {
      const gl = this.gl;

      const vertexShader = this.loadShader(gl.VERTEX_SHADER, vertexShaderSource);
      const fragmentShader = this.loadShader(gl.FRAGMENT_SHADER, fragmentShaderSource);

      const shaderProgram = gl.createProgram();
      gl.attachShader(shaderProgram, vertexShader);
      gl.attachShader(shaderProgram, fragmentShader);
      gl.linkProgram(shaderProgram);

      if(!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
        alert('Unable to initialize the shader program: ' + gl.getProgramInfoLog(shaderProgram));
        return null;
      }

      return shaderProgram;
    }

    attachTextures(framebuffer, textures) {
      const gl = this.gl;

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, framebuffer);

      const draw_buffers = [];

      for (let i = 0; i < textures.length; ++i) {
        gl.framebufferTexture2D(
          gl.DRAW_FRAMEBUFFER,
          gl['COLOR_ATTACHMENT' + i],
          gl.TEXTURE_2D,
          textures[i],
          0,
        );

        draw_buffers.push(gl['COLOR_ATTACHMENT' + i]);
      }

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);

      return draw_buffers;
    }

    initDefaultVertexBuffer() {
      const gl = this.gl;

      const vertex_array = new Float32Array([
        1, 1, 0,
        0, 1, 0,
        1, 0, 0,
        0, 0, 0,
      ]);

      const vertex_buffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, vertex_buffer);
      gl.bufferData(gl.ARRAY_BUFFER, vertex_array, gl.STATIC_DRAW, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, null);

      this.default_vertex_buffer = vertex_buffer;
    }

    useDefaultVertexBuffer(program) {
      const gl = this.gl;

      gl.bindBuffer(gl.ARRAY_BUFFER, this.default_vertex_buffer);

      const vertex_loc = gl.getAttribLocation(program, 'position');
      const numComponents = 3;
      const type = gl.FLOAT;
      const normalize = false;
      const stride = 0;
      const offset = 0;

      // The array buffer must be bound from earlier
      gl.vertexAttribPointer(
        vertex_loc,
        numComponents,
        type,
        normalize,
        stride,
        offset,
      );

      gl.enableVertexAttribArray(vertex_loc);

      gl.bindBuffer(gl.ARRAY_BUFFER, null);
    }

    getUintTextureArray(texture, width, height, array) {
      const gl = this.gl;

      const framebuffer = gl.createFramebuffer();

      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, framebuffer);
      gl.framebufferTexture2D(gl.READ_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
      gl.readBuffer(gl.COLOR_ATTACHMENT0);

      gl.readPixels(0, 0, width, height, gl.RGBA_INTEGER, gl.UNSIGNED_INT, array);

      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    }

    getFloatTextureArray(texture, width, height, array) {
      const gl = this.gl;

      const framebuffer = gl.createFramebuffer();

      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, framebuffer);
      gl.framebufferTexture2D(gl.READ_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
      gl.readBuffer(gl.COLOR_ATTACHMENT0);

      gl.readPixels(0, 0, width, height, gl.RGBA, gl.FLOAT, array);

      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    }

    runProgram(gl_helper, framebuffer, program, uniforms, locations, out_textures) {
      const gl = gl_helper.gl;

      gl.useProgram(program);
      gl_helper.setUniforms(uniforms, locations);

      const draw_buffers = gl_helper.attachTextures(framebuffer, out_textures);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, framebuffer);
      gl.drawBuffers(draw_buffers);

      gl_helper.useDefaultVertexBuffer(program);

      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);
    }

    runBigProgram(gl_helper, framebuffer, program, uniforms, locations, out_textures) {
      const gl = gl_helper.gl;
      const canvas = gl_helper.canvas;

      const old_width = canvas.width;
      const old_height = canvas.height;

      canvas.width = 2*old_width;
      canvas.height = 2*old_height;

      gl.viewport(0, 0, canvas.width, canvas.height);

      gl.useProgram(program);
      gl_helper.setUniforms(uniforms, locations);

      const draw_buffers = gl_helper.attachTextures(framebuffer, out_textures);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, framebuffer);
      gl.drawBuffers(draw_buffers);

      gl_helper.useDefaultVertexBuffer(program);

      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);

      canvas.width = old_width;
      canvas.height = old_height;

      gl.viewport(0, 0, canvas.width, canvas.height);
    }

    runSimulation(gl_helper, framebuffer, program, uniforms, locations, out_textures, cl_idx, clear) {
      const gl = gl_helper.gl;

      gl.useProgram(program);
      gl_helper.setUniforms(uniforms, locations, cl_idx);

      const draw_buffers = gl_helper.attachTextures(framebuffer, out_textures);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, framebuffer);
      gl.drawBuffers(draw_buffers);

      gl_helper.useDefaultVertexBuffer(program);

      if (clear) {
        gl.clearColor(0.0, 0.0, 0.0, 0.0);
        gl.clearDepth(1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      }

      gl.enable(gl.BLEND);
      gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE, gl.ONE, gl.ONE);
      gl.blendEquation(gl.FUNC_ADD);

      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

      // Reset defaults
      gl.blendFunc(gl.ONE, gl.ZERO);
      gl.disable(gl.BLEND);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);
    }

    runFinal(gl_helper, framebuffer, program, uniforms, locations, out_textures, cl_idx, simulation_length) {
      const gl = gl_helper.gl;
      const canvas = gl_helper.canvas;

      const old_width = canvas.width;
      const old_height = canvas.height;

      canvas.width = simulation_length;
      canvas.height = 1;

      gl.viewport(0, 0, canvas.width, canvas.height);

      gl.useProgram(program);
      gl_helper.setUniforms(uniforms, locations, cl_idx);

      const draw_buffers = gl_helper.attachTextures(framebuffer, out_textures);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, framebuffer);
      gl.drawBuffers(draw_buffers);

      gl_helper.useDefaultVertexBuffer(program);

      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);

      canvas.width = old_width;
      canvas.height = old_height;

      gl.viewport(0, 0, canvas.width, canvas.height);
    }

    setupDefault(shader_desc, context) {
      const gl = this.gl;

      const framebuffer = gl.createFramebuffer();
      const program = this.loadShaderProgram(shader_desc.vert, shader_desc.frag);
      gl.useProgram(program);

      const uniforms = shader_desc.uniforms;
      const locations = uniforms.map(u => gl.getUniformLocation(program, u[0]));

      const out_textures = shader_desc.out;

      // The bind is necessary to retain the correct context when the function is called
      const run = shader_desc.run.bind(context);

      return (...restArgs) => {
        run(this, framebuffer, program, uniforms, locations, out_textures, ...restArgs);
      };
    }

    setUniforms(uniforms, locations, ...uniformArgs) {
      const gl = this.gl;
      let texnum = 0;

      for (let i = 0; i < uniforms.length; ++i) {
        switch (uniforms[i][1]) {
          case '1f':
            gl.uniform1f(locations[i], uniforms[i][2](...uniformArgs));
            break;
          case '1i':
            gl.uniform1i(locations[i], uniforms[i][2](...uniformArgs));
            break;
          case '4fv_a':
            gl.uniform4fv(locations[i], ...uniforms[i][2](...uniformArgs));
            break;
          case 'tex':
            gl.activeTexture(gl['TEXTURE'+String(texnum)]);
            gl.bindTexture(gl.TEXTURE_2D, uniforms[i][2](...uniformArgs));
            gl.uniform1i(locations[i], texnum++);
            break;
          case 'mat4':
            gl.uniformMatrix4fv(locations[i], false, uniforms[i][2](...uniformArgs));
            break;
        }
      }
    }
  };
});
