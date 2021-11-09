import { join } from "path";

const resolve = (path: string) => join(__dirname, "./warthog/bin/", path);

export const warthog = resolve("./warthog");

export const roadhog = resolve("./roadhog");
