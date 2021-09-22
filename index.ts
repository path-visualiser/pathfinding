import { join } from "path";

const resolve = (path: string) => join(__dirname, path);

export const warthog = resolve("./warthog/bin/warthog");
