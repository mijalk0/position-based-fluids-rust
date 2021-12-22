#version 330
in vec3 position;
in vec3 normal;
in vec3 world_position;
in vec3 velocity;
in vec3 vel_color;

uniform mat4 view_matrix;
uniform mat4 proj_matrix;

out vec3 frag_position;
out vec3 frag_normal;
out vec3 frag_velocity;
out vec3 frag_vel_color;

void main() {
    vec3 final_position = position + world_position;
    frag_position =  final_position;
    frag_normal = normal;
    frag_velocity = velocity;
    frag_vel_color = vel_color;

    gl_Position = proj_matrix * view_matrix * vec4(final_position, 1.0);
}
