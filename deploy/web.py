"""
Entry point of microbeFunction_mcp — mounts all MCP servers under one FastAPI app.

Usage:
    cp default.conf.toml local.conf.toml
    uv run -m deploy.web
"""

import contextlib
import logging
from typing import List

import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from mcp.server.fastmcp import FastMCP

from tools.kegg.server import mcp as kegg_mcp
from tools.mgnify.server import mcp as mgnify_mcp

from deploy.config import conf

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

mcps: List[FastMCP] = [
    kegg_mcp,
    mgnify_mcp,
]


@contextlib.asynccontextmanager
async def lifespan(app: FastAPI):
    async with contextlib.AsyncExitStack() as stack:
        for mcp in mcps:
            await stack.enter_async_context(mcp.session_manager.run())
        yield


app = FastAPI(title="microbeFunction_mcp", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


for mcp in mcps:
    app.mount(f"/{mcp.name}/", mcp.streamable_http_app())


@app.get("/api/list_mcps")
def list_mcps():
    base_url = f"http://127.0.0.1:{conf['port']}"
    return {
        mcp.name: {
            "transport": "streamable_http",
            "url": f"{base_url}/{mcp.name}/mcp/",
        }
        for mcp in mcps
    }


if __name__ == "__main__":
    uvicorn.run(
        "deploy.web:app",
        host="0.0.0.0",
        workers=conf.get("workers", 1),
        port=conf["port"],
        reload=False,
    )
