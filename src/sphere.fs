#version 330
in vec3 frag_position;
in vec3 frag_normal;
in vec3 frag_velocity;
in vec3 frag_vel_color;

const float ambient_light = 0.4;
const vec3 light_dir = normalize(vec3(20.0, 30.0, 40.0));

out vec4 color;

void main() {
    float specular = max(dot(light_dir, frag_normal), ambient_light);
    vec3 scaled_color = specular * frag_vel_color;

    color = vec4(scaled_color, 1.0);
}
